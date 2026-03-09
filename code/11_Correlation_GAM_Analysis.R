# - create overview plot with CWSI, ecophysiological, and environmental
#   variables during the course of the experiment per treatment and tree species
# - calculate correlations between CWSI and environmental and ecophysiological
#   drought stress variables
# - create variable distribution overviews
# - create GAM Models to predict TWDmin classes (TWDmin = 0 and > 0) from
#   various variables
# - perform cross validation of GAM model
# - calculate evaluation matrices of GAM Models
# - calculate thresholds for 50% probability that TWDmin > 0 per tree species
#   and for various GAMs and visualize results per tree species
# - calculate cwsi means and sd betweeen TWDmin classes
# - create overview boxplots for CWSI, SWC and VPD values according to
#   TWDmin classes

library(readxl)
library(caret)
library(ggplot2)
library(glmmTMB)
library(gridExtra)
library(readr)
library(tidyr)
library(purrr)
library(corrplot)
library(lme4)
library(mgcv)
library(mgcViz)
library(broom.mixed) 
library(scales)
library(DHARMa)
library(car)
library(dplyr)
library(GGally)
library(visreg)
library(ggpubr)
library(patchwork)
library(dplyr)
library(scales)
library(lubridate)
library(patchwork)
library(zoo)

source("./Code/Functions.R")

irrig_CWSI <- read_csv("./data/results/Irrig_CWSI.csv")
# irrig_CWSI <- read_csv("./data/results/Irrig_CWSI_VPD_uncor.csv") # version with uncorrected VPD, airtemp, RH
# irrig_CWSI <- read_csv("./data/results/Irrig_CWSI_VPD_median_cor.csv") # version with Tmedian corrected VPD, airtemp, RH

all_data10 <- read_csv("./data/results/Env_Eco_CWSI_combined_10min.csv")
# all_data10 <- read_csv("./data/results/Env_Eco_CWSI_combined_10min_VPD_uncor.csv")
# all_data10 <- read_csv("./data/results/Env_Eco_CWSI_combined_10min_VPD_median_cor.csv")

# create overview plot with CWSI, ecophysiological, and environmental
# variables during the course of the experiment per treatment and tree species

# calculate daily values
all_data_stats <- all_data10 %>%
  group_by(Date, Treatment, Tree.Species) %>%
  summarise(mean_Min_rel_TWD  = mean(Min_rel_TWD, na.rm = TRUE),
            sd_Min_rel_TWD    = sd(Min_rel_TWD, na.rm = TRUE),
            
            mean_Max_rel_TWD  = mean(Max_rel_TWD, na.rm = TRUE),
            sd_Max_rel_TWD    = sd(Max_rel_TWD, na.rm = TRUE),
            
            mean_rel_TWD = mean(rel_TWD_log, na.rm = TRUE),
            sd_rel_TWD   = sd(rel_TWD_log, na.rm = TRUE),
            
            mean_SWC = mean(VWC_mean, na.rm = TRUE),
            sd_SWC   = sd(VWC_mean, na.rm = TRUE),
            
            mean_WP_midday = mean(WP_midday, na.rm = TRUE),
            sd_WP_midday   = sd(WP_midday, na.rm = TRUE),
            
            mean_WP_predawn = mean(WP_predawn, na.rm = TRUE),
            sd_WP_predawn   = sd(WP_predawn, na.rm = TRUE),
            
            mean_CWSI = mean(CWSI, na.rm = TRUE),
            sd_CWSI   = sd(CWSI, na.rm = TRUE),
            
            mean_VPD_kPa = mean(VPD_kPa, na.rm = TRUE),
            sd_VPD = sd(VPD_kPa, na.rm = TRUE),
            min_VPD_kPa = min(VPD_kPa, na.rm = TRUE),
            max_VPD_kPa = max(VPD_kPa, na.rm = TRUE)
            )

all_data_filled <- all_data_stats %>%
  group_by(Tree.Species, Treatment) %>%
  complete(Date = seq.Date(min(Date), max(Date), by = "day")) %>%
  arrange(Date) %>%
  mutate(
    mean_CWSI = na.approx(mean_CWSI, Date, na.rm = FALSE),
    sd_CWSI   = na.approx(sd_CWSI, Date, na.rm = FALSE),
    mean_WP_midday = na.approx(mean_WP_midday, Date, na.rm = FALSE),
    mean_WP_predawn = na.approx(mean_WP_predawn, Date, na.rm = FALSE),
  ) %>%
  ungroup()


colorvalues <- c(
  "Fagus sylvatica" = "#298c8c",
  "Pseudotsuga menziesii" = "orange",
  "Quercus robur"         = "darkmagenta"
)

CWSI_plot <- ggplot() +
  geom_line(data = all_data_filled,
            aes(x = Date,
                y = mean_CWSI,
                color = Tree.Species,
                linetype = Treatment,
                group = interaction(Tree.Species, Treatment)),
            linewidth = 0.5) +
  geom_point(data = all_data_stats,
             aes(x = Date, 
                 y = mean_CWSI, 
                 color = Tree.Species,
                 shape = Treatment,
                 group = interaction(Tree.Species, Treatment)),
             size = 2) +  
  scale_color_manual(values = colorvalues) +
  scale_linetype_manual(values = c(
    "Watered"   = "dotdash",
    "Droughted" = "solid"
  )) +
  scale_shape_manual(values = c(
    "Watered"   = 16,   
    "Droughted" = 17    
  )) +
  labs(y = "CWSI",
       x = "Date") +
  xlim(as.Date("2023-06-25"), as.Date("2023-08-20")) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),  
    axis.text.x  = element_blank(),  
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )

TWDmin_plot <- ggplot(all_data_stats) +
  geom_line(aes(x = Date, 
                y = mean_Min_rel_TWD, 
                color = Tree.Species,
                linetype = Treatment,       
                group = interaction(Tree.Species, Treatment)),
            linewidth = 0.8) +
  
  scale_color_manual(values = colorvalues) +
  
  scale_linetype_manual(values = c(
    "Watered" = "dotdash",
    "Droughted" = "solid"
  )) +
  labs(y = expression(TWD[min]))+
  xlim(as.Date("2023-06-25"), as.Date("2023-08-20")) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),  
    axis.text.x  = element_blank(),  
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )


TWDmax_plot <- ggplot(all_data_stats) +
  geom_line(aes(x = Date, 
                y = mean_Max_rel_TWD, 
                color = Tree.Species,
                linetype = Treatment,       
                group = interaction(Tree.Species, Treatment)),
            linewidth = 0.8) +
  
  scale_color_manual(values = colorvalues,
                     name = "Tree species") +
  
  scale_linetype_manual(values = c(
    "Watered" = "dotdash",
    "Droughted" = "solid"
  )) +
  labs(y = expression(TWD[max]))+
  xlim(as.Date("2023-06-25"), as.Date("2023-08-20")) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),  
    axis.text.x  = element_blank(),  
    axis.ticks.x = element_blank()   
  )

TWDrel_plot <- ggplot(all_data_stats) +
  geom_line(aes(x = Date, 
                y = mean_rel_TWD, 
                color = Tree.Species,
                linetype = Treatment,       
                group = interaction(Tree.Species, Treatment)),
            linewidth = 0.8) +
  
  scale_color_manual(values = colorvalues) +
  
  scale_linetype_manual(values = c(
    "Watered" = "dotdash",
    "Droughted" = "solid"
  )) +
  labs(y = expression(TWD[rel]))+
  xlim(as.Date("2023-06-25"), as.Date("2023-08-20")) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),  
    axis.text.x  = element_blank(),  
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )

PDWP_plot <- ggplot() +
  geom_line(data = all_data_filled,
            aes(x = Date,
                y = mean_WP_predawn,
                color = Tree.Species,
                linetype = Treatment,
                group = interaction(Tree.Species, Treatment)),
            linewidth = 0.5) +
  geom_point(data = all_data_stats,
             aes(x = Date, 
                 y = mean_WP_predawn, 
                 color = Tree.Species,
                 shape = Treatment,
                 group = interaction(Tree.Species, Treatment)),
             size = 2) +  
  scale_color_manual(values = colorvalues) +
  scale_linetype_manual(values = c(
    "Watered"   = "dotdash",
    "Droughted" = "solid"
  )) +
  scale_shape_manual(values = c(
    "Watered"   = 16,   
    "Droughted" = 17    
  )) +
  labs(y = bquote(psi["predawn"] ~ "[MPa]"),
       x = "Date") +
  xlim(as.Date("2023-06-25"), as.Date("2023-08-20")) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),  
    axis.text.x  = element_blank(),  
    axis.ticks.x = element_blank()   
  )

MDWP_plot <- ggplot() +
  geom_line(data = all_data_filled,
            aes(x = Date,
                y = mean_WP_midday,
                color = Tree.Species,
                linetype = Treatment,
                group = interaction(Tree.Species, Treatment)),
            linewidth = 0.5) +
  geom_point(data = all_data_stats,
             aes(x = Date, 
                 y = mean_WP_midday, 
                 color = Tree.Species,
                 shape = Treatment,
                 group = interaction(Tree.Species, Treatment)),
             size = 2) +  
  scale_color_manual(values = colorvalues) +
  scale_linetype_manual(values = c(
    "Watered"   = "dotdash",
    "Droughted" = "solid"
  )) +
  scale_shape_manual(values = c(
    "Watered"   = 16,   
    "Droughted" = 17    
  )) +
  labs(y = bquote(psi["midday"] ~ "[MPa]"),
       x = "Date") +
  xlim(as.Date("2023-06-25"), as.Date("2023-08-20")) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),  
    axis.text.x  = element_blank(),  
    axis.ticks.x = element_blank() ,
    legend.position = "none"
  )
  

SWC_plot <- ggplot(all_data_stats) +
  geom_line(aes(x = Date, 
                y = mean_SWC, 
                color = Tree.Species,
                linetype = Treatment,       
                group = interaction(Tree.Species, Treatment)),
            linewidth = 0.8) +
  scale_color_manual(values = colorvalues) +
  scale_linetype_manual(values = c(
    "Watered" = "dotdash",
    "Droughted" = "solid"
  )) +
  labs(y = "SWC [%]")+
  xlim(as.Date("2023-06-25"), as.Date("2023-08-20")) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),  
    axis.text.x  = element_blank(),  
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )

VPD_plot <- ggplot(all_data_stats) +
  geom_ribbon(aes(x = Date, ymin = min_VPD_kPa, ymax = max_VPD_kPa, 
                  fill = "VPD range")) +
  geom_line(aes(x = Date, y = mean_VPD_kPa, color = "Mean VPD"),
            linewidth = 0.8) +
  scale_fill_manual(name = "", values = c("VPD range" = "#a6cee3")) +
  scale_color_manual(name = "", values = c("Mean VPD" = "#1f78b4")) +
  labs(y = "VPD [kPa]", x = "2023") +
  xlim(as.Date("2023-06-25"), as.Date("2023-08-20")) +
  theme_bw() +
  theme()



# combine all plots
(CWSI_plot / TWDmin_plot / TWDmax_plot / TWDrel_plot / PDWP_plot / MDWP_plot /
    SWC_plot / VPD_plot) 


ggsave("./graphics/results/Dataoverview_Treatment_Species.png",
       width = 21,      
       height = 29,     
       units = "cm",    
       dpi = 300)

# ggsave("./graphics/results/Dataoverview_Treatment_Species_VPD_uncor.png",
#        width = 21,
#        height = 29,
#        units = "cm",
#        dpi = 300)
# 
# ggsave("./graphics/results/Dataoverview_Treatment_Species_VPD_median_cor.png",
#        width = 21,      
#        height = 29,     
#        units = "cm",    
#        dpi = 300)

ggsave("./graphics/results/Dataoverview_Treatment_Species2.png",
       width = 25,      
       height = 29,     
       units = "cm",    
       dpi = 300)

# Correlation analysis
vars <- irrig_CWSI %>%
  dplyr::select(Mean_VPD_kPa, Min_VPD_kPa, Max_VPD_kPa, VPD_kPa, RH, Air_temp,
         SWC,
         Min_rel_TWD,         
         Max_rel_TWD, WP_midday, WP_predawn, rel_TWD_log) %>%
  names()

cor_spearman <- function(x, y) {
  df <- data.frame(x = x, y = y) %>% na.omit()
  if (nrow(df) < 3) {
    return(tibble(correlation = NA_real_, p_value = NA_real_))
  }
  test <- cor.test(df$x, df$y, method = "spearman", exact = FALSE)
  tibble(correlation = unname(test$estimate), p_value = test$p.value)
}

# Using all variables in the correlation analysis, perform by tree species
cor_per_species <- irrig_CWSI %>%
  filter(!is.na(CWSI)) %>%
  group_by(Tree.Species) %>%
  group_map(~ {
    # .x is the subset of the dataframe per tree species
    map_dfr(vars, function(var) {
      res <- cor_spearman(.x$CWSI, .x[[var]])
      tibble(
        Tree.Species = unique(.x$Tree.Species),
        variable = var,
        correlation = res$correlation,
        p_value = res$p_value
      )
    })
  }, .keep = TRUE) %>%
  bind_rows() %>%
  mutate(
    p_value = ifelse(p_value < 0.001, "< 0.001", sprintf("%.3f", p_value)),
    correlation = round(correlation, 2)
  )

# calculate correlation across all tree species
cor_overall <- map_dfr(vars, function(var) {
  cor_spearman(irrig_CWSI$CWSI, irrig_CWSI[[var]]) %>%
    mutate(variable = var)
}) %>%
  mutate(
    p_value = ifelse(p_value < 0.001, "< 0.001", sprintf("%.3f", p_value)),
    correlation = round(correlation, 2),
    Tree.Species = "Overall"
  )

correlation_species_df <- rbind(cor_overall, cor_per_species)

write.csv(correlation_species_df, 
          "./data/results/statistical_analysis/Correlations_BDE.csv",
          row.names = FALSE)

# write.csv(correlation_species_df, 
#           "./data/results/statistical_analysis_VPD_uncor/Correlations_BDE_VPD_uncor.csv",
#           row.names = FALSE)
# 
# write.csv(correlation_species_df,
#           "./data/results/statistical_analysis_VPD_median_cor/Correlations_BDE_VPD_median_cor.csv",
#           row.names = FALSE)

formatted_df <- correlation_species_df %>%
  mutate(
    correlation = sprintf("%.2f", correlation),
    cor_p = paste0(correlation, " (", p_value, ")")
  ) %>%
  dplyr::select(variable, Tree.Species, cor_p)

correlation_species_df_wide <- formatted_df %>%
  pivot_wider(
    names_from = Tree.Species,
    values_from = cor_p
  )

write.csv(correlation_species_df_wide, 
          "./data/results/statistical_analysis/Correlations_BDE_formatted.csv",
          row.names = FALSE)

# write.csv(correlation_species_df_wide, 
#           "./data/results/statistical_analysis_VPD_uncor/Correlations_BDE_formatted_VPD_uncor.csv",
#           row.names = FALSE)
# 
# write.csv(correlation_species_df_wide,
#           "./data/results/statistical_analysis_VPD_median_cor/Correlations_BDE_formatted_VPD_median_cor.csv",
#           row.names = FALSE)

corplot_data <- irrig_CWSI %>%
  filter(!is.na(CWSI)) %>%
  dplyr::select(CWSI, all_of(vars))

mat<- cor(corplot_data, use = "pairwise.complete.obs", method = "spearman")
corrplot(mat, method = "color")

# Calculate p-values for correlation
cor.mtest <- function(data, conf.level = 0.95) {
  data <- as.matrix(data)
  n <- ncol(data)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      tmp <- cor.test(data[,i], data[,j], method = "spearman")
      p.mat[i,j] <- tmp$p.value
      p.mat[j,i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(data)
  return(p.mat)
}

p.mat <- cor.mtest(corplot_data)

# Create corrplot with significance
corrplot(mat, method = "color",
         type = "upper",
         order = "hclust",
         addCoef.col = "black",
         tl.col = "black", tl.srt = 45,
         p.mat = p.mat, sig.level = 0.01, insig = "blank",
         diag = FALSE)


# correlation matrix (spearman)
cor_matrix <- cor(corplot_data, use = "pairwise.complete.obs", 
                  method = "spearman")

# Extract correlation values (besides correlation between CWSI and itself)
cwsi_corr <- cor_matrix["CWSI", ]
cwsi_corr <- cwsi_corr[names(cwsi_corr) != "CWSI"]

cwsi_corr_df <- data.frame(
  Variable = names(cwsi_corr),
  Correlation = cwsi_corr
)

# Calculate p-values for correlation
cor_pvals <- sapply(names(cwsi_corr), function(var) {
  cor.test(corplot_data[[var]], corplot_data$CWSI, method = "spearman")$p.value
})
cwsi_corr_df$p_value <- cor_pvals

# add significance levels
cwsi_corr_df$Significant <- ifelse(cwsi_corr_df$p_value < 0.001, "***",
                                   ifelse(cwsi_corr_df$p_value < 0.01, "**",
                                       ifelse(cwsi_corr_df$p_value < 0.05, "*",
                                                 "")))
# Plot correlation across all tree species
ggplot(
  cwsi_corr_df %>%
    mutate(
      Variable = dplyr::recode(
        Variable,
        "Mean_VPD_kPa" = "VPD[mean]",
        "Min_VPD_kPa"  = "VPD[min]",
        "Max_VPD_kPa"  = "VPD[max]",
        "VPD_kPa"      = "VPD",
        "RH"           = "RH",
        "Air_temp"     = "T[air]",
        "SWC"          = "SWC",
        "Min_rel_TWD"  = "TWD[min]",
        "Max_rel_TWD"  = "TWD[max]",
        "WP_midday"    = "psi[midday]",
        "WP_predawn"   = "psi[predawn]",
        "rel_TWD_log"  = "TWD[rel]"
      )
    ),
  aes(
    y = reorder(Variable, -abs(Correlation)),
    x = Correlation,
    fill = Correlation
  )
) +
  geom_bar(stat = "identity") +
  geom_text(
    aes(label = Significant),
    size = 5,
    hjust = 0.5,
    vjust = ifelse(cwsi_corr_df$Correlation >= 0, -0.5, 1.5)
  ) +
  coord_flip() +
  scale_y_discrete(labels = label_parse()) +   
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Correlation of variables with CWSI",
    x = "Spearman’s rho (ρ)",
    y = "Variable",
    fill = "ρ"
  ) +
  xlim(-0.6, 0.6) +
  theme(
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text.x  = element_text(size = 15),
    axis.text.y  = element_text(size = 15),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text  = element_text(size = 15)
  )

ggsave("./graphics/results/Corr_Plot_CWSI_all_trees_BED.png",  
       width = 14, height = 6, dpi = 300)

# Plot correlations with CWSI divided into ecophysiological and 
# environmental variables (Figure 3)
ggplot(
  cwsi_corr_df %>%
    mutate(
      Variable_cat = case_when(
        Variable %in% c("Min_rel_TWD", "Max_rel_TWD", "WP_midday", "WP_predawn",
                        "rel_TWD_log") ~ "Ecophysiological variables",
        TRUE ~ "Environmental variables"),
      Variable = dplyr::recode(
        Variable,
        "Mean_VPD_kPa" = "VPD[mean]",
        "Min_VPD_kPa"  = "VPD[min]",
        "Max_VPD_kPa"  = "VPD[max]",
        "VPD_kPa"      = "VPD",
        "RH"           = "RH",
        "Air_temp"     = "T[air]",
        "SWC"          = "SWC",
        "Min_rel_TWD"  = "TWD[min]",
        "Max_rel_TWD"  = "TWD[max]",
        "WP_midday"    = "psi[midday]",
        "WP_predawn"   = "psi[predawn]",
        "rel_TWD_log"  = "TWD[rel]"
      )
    ),
  aes(
    y = reorder(Variable, -abs(Correlation)),
    x = Correlation,
    fill = Correlation
  )
) +
  geom_bar(stat = "identity") +
  geom_text(
    aes(label = Significant),
    size = 5,
    hjust = 0.5,
    vjust = ifelse(cwsi_corr_df$Correlation >= 0, -0.5, 1.5)
  ) +
  coord_flip() +
  scale_y_discrete(labels = label_parse()) +   
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)+
  theme_bw(base_size = 14) +
  labs(
    x = "Spearman’s rho (ρ)",
    y = "",
    fill = "ρ"
  ) +
  xlim(-0.6, 0.6) +
  facet_grid(. ~ Variable_cat, scale = "free") +
  theme(
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text.x  = element_text(size = 15),
    axis.text.y  = element_text(size = 15),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text  = element_text(size = 15),
    strip.text = element_text(size = 20, face = "bold")
  )

ggsave("./graphics/results/Corr_Plot_CWSI_all_trees_BED_env_eco2.png",  
       width = 14, height = 6, dpi = 300)
# ggsave("./graphics/results/Corr_Plot_CWSI_all_trees_BED_env_eco2_VPD_uncor.png",  
#        width = 14, height = 6, dpi = 300)
# ggsave("./graphics/results/Corr_Plot_CWSI_all_trees_BED_env_eco2_VPD_median_cor.png",
#        width = 14, height = 6, dpi = 300)


# Analyse dendrometer trees
data_B <- irrig_CWSI %>%
  filter(Dendro_info == "DendroTree" & !is.na(CWSI), !is.na(Min_rel_TWD), 
         !is.na(SWC), !is.na(VPD_kPa))

n_WP_midday <- data_B %>%
  filter(!is.na(WP_midday)) %>%
  distinct(Date)%>%
  nrow()

n_WP_predawn <- data_B %>%
  filter(!is.na(WP_predawn)) %>%
  distinct(Date)%>%
  nrow()

# plot histogramms per tree species for the predictor variables, MirTWD 
# and relTWD
vars_to_plot <- c("SWC", "VPD_kPa", "Min_rel_TWD", "CWSI", "rel_TWD_log")

tree_species <- unique(data_B$Tree.Species)

for(var in vars_to_plot){
  p <- ggplot(data_B, aes(x = .data[[var]])) +
    geom_histogram(bins = 60, fill = "skyblue", color = "black") +
    theme_minimal() +
    facet_wrap(~ Tree.Species, scales = "free") +  
    ggtitle(paste("Histogram of", var, "by Tree Species")) +
    xlab(var) +
    ylab("Count") +
    theme(plot.title = element_text(hjust = 0.5))

  filename <- paste0("./graphics/histogramms_DendroTrees/", 
                     var, "_by_TreeSpecies.png")
  ggsave(filename, plot = p, width = 10, height = 6, dpi = 300)
}


# correlation plots
for(tree_spec in tree_species){
  data_sub <- data_B %>%
    filter(Tree.Species == tree_spec) %>%
    dplyr::select(Min_rel_TWD, CWSI, SWC, VPD_kPa) %>%  
    rename(
      TWDmin = Min_rel_TWD,
      # TWDrel = rel_TWD_log,
      CWSI = CWSI,
      SWC = SWC,
      VPD = VPD_kPa
    ) %>%
    na.omit()  
  
  # ggpairs Correlation plot
  p <- ggpairs(data_sub,
               upper = list(continuous = wrap("cor", size = 4,
                                              method = "spearman")),
               lower = list(continuous = "points"),
               diag = list(continuous = "barDiag"),
               ) +
    theme_bw()
  print(p)
  
    ggsave(filename = paste0("./graphics/results/correlations/corrplot_",
                           tree_spec, ".png"),
         plot = p, width = 8, height = 8)
    
}


# Overview with data amount (0/non for MirTWD)
summary_table_species <- data_B %>%
  group_by(Tree.Species) %>%
  summarise(
    SWC = sum(!is.na(SWC)),
    VPD_kPa = sum(!is.na(VPD_kPa)),
    Min_rel_TWD_Total = sum(!is.na(Min_rel_TWD)),
    CWSI = sum(!is.na(CWSI)),
    rel_TWD_log = sum(!is.na(rel_TWD_log)),
    Min_rel_TWD_Zero = sum(Min_rel_TWD == 0, na.rm = TRUE),
    Min_rel_TWD_NonZero = sum(Min_rel_TWD != 0, na.rm = TRUE)
  )


write.csv(summary_table_species,
          "./data/results/statistical_analysis/Value_Overview_Species_Variable_BDE.csv")

# Prediction of TWDmin classes (TWDmin =0 and > 0)
set.seed(123)

data_B <- irrig_CWSI %>%
  filter(Dendro_info == "DendroTree" & !is.na(CWSI),
         !is.na(Min_rel_TWD), !is.na(SWC), !is.na(VPD_kPa)) %>%
  dplyr::select(Tree.Species, CWSI, Min_rel_TWD, SWC, VPD_kPa) %>%
  unique() %>%
  mutate(class = factor(ifelse(Min_rel_TWD == 0, "zero", "nonzero"),
                        levels = c("zero", "nonzero")),
         Tree.Species = factor(Tree.Species))

# check normality
par(mfrow = c(2,2))
for(i in c(2:5)){
  qqnorm(data_B[[i]], main = names(data_B)[i])
  qqline(data_B[[i]])
}
par(mfrow = c(1,1))

# correlations across all tree species
p <- ggpairs(data_B %>%
          dplyr::select(-c("Tree.Species", "class")) %>%
          rename(
            TWDmin = Min_rel_TWD,
            
            CWSI = CWSI,
            SWC = SWC,
            VPD = VPD_kPa
          ) %>%
          na.omit() ,
        upper = list(continuous = wrap("cor", size = 4,
                                       method = "spearman")),
        lower = list(continuous = "points"),
        diag = list(continuous = "barDiag")) +
  theme_bw()

print(p)

ggsave(filename = "./graphics/results/correlations/corrplot_BED_overall.png",
       width = 10, height = 10)


# GAM Analysis
# Create predictor sets
predictor_sets <- list(
  CWSI_only       = c("CWSI"),
  CWSI_SWC        = c("CWSI", "SWC"),
  CWSI_VPD        = c("CWSI", "VPD_kPa"),
  CWSI_SWC_VPD    = c("CWSI", "SWC", "VPD_kPa"),
  SWC_only        = c("SWC"),
  VPD_only        = c("VPD_kPa"),
  SWC_VPD         = c("VPD_kPa", "SWC")
)


# Crossvalidation per tree species with random folds (OOF-Evaluation)

set.seed(123)

# Create stratified folds
k <- 5
folds <- createFolds(data_B$class, k = k, list = TRUE)

# Create folder for confusion matrix plots
dir.create("./graphics/Confusion_Matrix/plots_cm_bal", 
           showWarnings = FALSE, recursive = TRUE)
dir.create("./graphics/Confusion_Matrix/plots_cm_mean_bal", 
           showWarnings = FALSE, recursive = TRUE)

# initialise results lists and dataframes
cv_results_all_bal <- list()
cm_mean_all_bal <- list()
oof_results_all_bal <- data.frame()


# Loop over all predictor sets and perform cross validation
for (pred_set_name in names(predictor_sets)) {
  
  features <- predictor_sets[[pred_set_name]]
  cat("\n===== Predictor-Set:", pred_set_name, "=====\n")
  
  metrics_gam <- NULL
  cm_folds_list <- list()
  obs_pred_all <- data.frame()
  
  # Loop per fold
  for (fold_idx in seq_along(folds)) {
    
    test_idx  <- folds[[fold_idx]]
    train_idx <- setdiff(seq_len(nrow(data_B)), test_idx)
    
    train <- data_B[train_idx, ]
    test  <- data_B[test_idx, ]
    
    cat(paste0("Fold ", fold_idx, ":\n"))
    cat("  Train counts:\n")
    print(table(train$class))
    cat("  Test counts:\n")
    print(table(test$class))
    cat("\n")
    
    # Apply GAM
    gam_formula <- as.formula(
      paste("class ~ s(", paste(features, collapse = ") + s("),
            ") + Tree.Species")
    )
    
    gam_model <- gam(gam_formula, family = binomial, data = train)
    
    # Predict (probabilities)
    gam_probs <- predict(gam_model, newdata = test, type = "response")
    
    # assign binary classes (TWdmin = 0 and > 0)
    gam_pred <- factor(ifelse(gam_probs >= 0.5, "nonzero", "zero"),
                       levels = levels(data_B$class))
    
    # Calculate performance metrics per fold
    metrics_fold <- extract_metrics(gam_pred, test$class)
    metrics_gam <- rbind(metrics_gam, metrics_fold)
    
    # confusion matrix per fold
    cm_folds_list[[fold_idx]] <- table(
      factor(test$class, levels = levels(data_B$class)),
      gam_pred
    )
    
    # Collect OOF-results
    obs_pred_all <- rbind(
      obs_pred_all,
      data.frame(Observed = test$class, Predicted = gam_pred, Model = "GAM")
    )
  }

  # Calculate mean cross-validation metrics across all folds
  metrics_gam_summary <- data.frame(
    Metric = names(colMeans(metrics_gam, na.rm = TRUE)),
    Mean   = colMeans(metrics_gam, na.rm = TRUE),
    SD     = apply(metrics_gam, 2, sd, na.rm = TRUE)
  )
  
  cv_results_all_bal[[pred_set_name]] <- list(GAM = metrics_gam_summary)
  
  # Calculate OOF-metrics 
  oof_results <- obs_pred_all %>%
    group_by(Model) %>%
    group_modify(~ {
      metrics <- extract_metrics(.x$Predicted, .x$Observed)
      as.data.frame(as.list(metrics))
    }) %>%
    ungroup()
  
  oof_results$Predictor_Set <- pred_set_name
  oof_results$Tree.Species  <- "all"
  oof_results$dist0N0       <- "equal_per_fold"
  
  oof_results_all_bal <- bind_rows(oof_results_all_bal, oof_results)
  
  # Plot OOF-confusionsmatrix 
  obs_classes <- sort(unique(obs_pred_all$Observed))
  pred_classes <- sort(unique(obs_pred_all$Predicted))
  
  all_combinations <- expand.grid(
    Model = unique(obs_pred_all$Model),
    Observed = obs_classes,
    Predicted = pred_classes
  )
  
  cm_data <- obs_pred_all %>%
    group_by(Model, Observed, Predicted) %>%
    summarise(n = n(), .groups = "drop") %>%
    right_join(all_combinations, by = c("Model", "Observed", "Predicted")) %>%
    mutate(n = replace_na(n, 0)) %>%
    group_by(Model, Observed) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup()
  
  color_map <- expand.grid(Observed = obs_classes, Predicted = pred_classes) %>%
    mutate(
      fill_color = case_when(
        Observed == obs_classes[1] & Predicted == pred_classes[1] ~ alpha("blue4", 0.7),
        Observed == obs_classes[1] & Predicted == pred_classes[2] ~ alpha("lightblue", 0.7),
        Observed == obs_classes[2] & Predicted == pred_classes[1] ~ alpha("pink", 0.7),
        Observed == obs_classes[2] & Predicted == pred_classes[2] ~ alpha("red4", 0.7)
      )
    )
  
  cm_data <- left_join(cm_data, color_map, by = c("Observed", "Predicted"))
  
  cm_plot <- ggplot(cm_data, aes(x = Observed, y = Predicted, fill = fill_color)) +
    geom_tile(color = "black") +
    geom_text(aes(label = paste0(n, " (", sprintf("%.1f%%", 100 * prop), ")")),
              fontface = "bold") +
    scale_fill_identity() +
    facet_wrap(~Model) +
    labs(
      title = paste("Out-of-Fold Confusion Matrix (%) -", pred_set_name),
      x = "Observed", y = "Predicted"
    ) +
    theme_minimal(base_size = 14)
  
  ggsave(
    filename = paste0("./graphics/Confusion_Matrix/plots_cm_bal/cm_oof_",
                      pred_set_name, ".png"),
    plot = cm_plot, width = 10, height = 6, dpi = 300
  )
  
  
  # Plot confusion matrix mean over folds
  all_levels <- levels(data_B$class)
  cm_sum <- matrix(0, nrow = length(all_levels), ncol = length(all_levels),
                   dimnames = list(Observed = all_levels, 
                                   Predicted = all_levels))
  for (cm in cm_folds_list) cm_sum <- cm_sum + as.matrix(cm)
  cm_mean <- cm_sum / length(cm_folds_list)
  cm_mean_df <- as.data.frame(as.table(cm_mean))
  cm_mean_df$Model <- "GAM"
  
  cm_mean_df <- cm_mean_df %>%
    group_by(Model, Observed) %>%
    mutate(prop = Freq / sum(Freq)) %>%
    ungroup() %>%
    left_join(color_map, by = c("Observed", "Predicted"))
  
  cm_plot_mean <- ggplot(cm_mean_df, aes(x = Observed, y = Predicted, 
                                         fill = fill_color)) +
    geom_tile(color = "black") +
    geom_text(aes(label = paste0(sprintf("%.1f", Freq), 
                                 " (", sprintf("%.1f%%", 100 * prop), ")")),
              fontface = "bold") +
    scale_fill_identity() +
    facet_wrap(~Model) +
    labs(
      title = paste("Average Confusion Matrix over folds -", pred_set_name),
      x = "Observed", y = "Predicted"
    ) +
    theme_minimal(base_size = 14)
  
  ggsave(
    filename = paste0("./graphics/Confusion_Matrix/plots_cm_mean_bal/cm_mean_",
                      pred_set_name, ".png"),
    plot = cm_plot_mean, width = 10, height = 6, dpi = 300
  )

  cm_mean_all_bal[[pred_set_name]] <- list(GAM = cm_mean_df)
}

# Collect results in long format
cv_summary_df_long_bal <- lapply(cv_results_all_bal, function(mod_list) {
  bind_rows(lapply(names(mod_list), function(mod_name) {
    df <- mod_list[[mod_name]]
    df$Model <- mod_name
    df
  }), .id = "Predictor_Set")
}) %>%
  bind_rows(.id = "Predictor_Set")

write.csv(cv_summary_df_long_bal, 
          "./data/results/statistical_analysis/Mean+SD_CV_results_GAM.csv")

# write.csv(cv_summary_df_long_bal,
#           "./data/results/statistical_analysis_VPD_uncor/Mean+SD_CV_results_GAM_VPD_uncor.csv")
# 
# write.csv(cv_summary_df_long_bal,
#           "./data/results/statistical_analysis_VPD_median_cor/Mean+SD_CV_results_GAM_VPD_median_cor.csv")

write.csv(oof_results_all_bal, 
          "./data/results/statistical_analysis/OOF_CV_results_GAM.csv")

# write.csv(oof_results_all_bal,
#           "./data/results/statistical_analysis_VPD_uncor/OOF_CV_results_GAM_VPD_uncor.csv")
# 
# write.csv(oof_results_all_bal,
#           "./data/results/statistical_analysis_VPD_median_cor/OOF_CV_results_GAM_VPD_median_cor.csv")


# Apply GAMs on whole dataset and diagnose with DHARMa 

# Create folders for DHARMa plots 
dir.create("./graphics/DHARMa_Plots", showWarnings = FALSE, recursive = TRUE)
dir.create("./data/results/statistical_analysis", showWarnings = FALSE, 
           recursive = TRUE)

# dir.create("./graphics/DHARMa_Plots_VPD_uncor", showWarnings = FALSE, recursive = TRUE)
# dir.create("./data/results/statistical_analysis_VPD_uncor", showWarnings = FALSE, 
#            recursive = TRUE)
# 
# 
# dir.create("./graphics/DHARMa_Plots_VPD_median_cor", showWarnings = FALSE, recursive = TRUE)
# dir.create("./data/results/statistical_analysis_VPD_median_cor", showWarnings = FALSE, 
#            recursive = TRUE)

# Create dataframes for results
dharma_summary <- data.frame()
model_summary <- data.frame()
gam_list <- list()

# Loop over all predictor-sets
for (pred_set_name in names(predictor_sets)) {
  
  features <- predictor_sets[[pred_set_name]]
  cat("\n===== GAM DHARMa Diagnostics:", pred_set_name, "=====\n")
  
  # Define GAM
  gam_formula <- as.formula(
    paste("class ~ s(", paste(features, collapse = ") + s("), 
          ") + Tree.Species")
  )
  
  # Fit model on whole dataframe
  gam_model <- gam(gam_formula, family = binomial, data = data_B)
  
  gam_list[[pred_set_name]] <- gam_model
  
  # DHARMa-Residualdiagnosis
  sim_res <- simulateResiduals(fittedModel = gam_model, plot = FALSE)
  
  # Combined DHARMa diagnosis plot
  png(filename = paste0("./graphics/DHARMa_Plots/DHARMa_Disp_ZI_",
                        pred_set_name, ".png"),
      width = 2400, height = 1600, res = 200)

  # png(filename = paste0("./graphics/DHARMa_Plots_VPD_uncor/DHARMa_Disp_ZI_",
  #                       pred_set_name, ".png"),
  #     width = 2400, height = 1600, res = 200)
  # 
  # png(filename = paste0("./graphics/DHARMa_Plots_VPD_median_cor/DHARMa_Disp_ZI_",
  #                       pred_set_name, ".png"),
  #     width = 2400, height = 1600, res = 200)
  
  par(mfrow = c(1, 2), mar = c(4,4,3,1)) 
  
  # Test for dispersion
  testDispersion(sim_res, plot = TRUE)
  
  # Test for Zero-Inflation
  testZeroInflation(sim_res, plot = TRUE)
  
  dev.off()
  
  png(filename = paste0("./graphics/DHARMa_Plots/DHARMa_plot_",
                        pred_set_name, ".png"),
      width = 2400, height = 1600, res = 200)

  # # DHARMa-Standardplots (Residuals vs Predicted + QQ)
  plot(sim_res)

  dev.off()
  
  ks <- testUniformity(sim_res)
  disp <- testDispersion(sim_res)
  zero <- testZeroInflation(sim_res)
  
  # Save DHARMa results
  dharma_summary <- rbind(
    dharma_summary,
    data.frame(
      Predictor_Set = pred_set_name,
      KS_Statistic  = ks$statistic,
      KS_pvalue     = ks$p.value,
      Dispersion_Statistic = disp$statistic,
      Dispersion_pvalue    = disp$p.value,
      ZeroInflation_Statistic = zero$statistic,
      ZeroInflation_pvalue    = zero$p.value
    )
  )
  
  # Plot effects
  png(filename = paste0("./graphics/DHARMa_Plots/Effects_",
                        pred_set_name, ".png"),
      width = 1800, height = 1000, res = 200)
  par(mfrow = c(ceiling(length(features)/2), 2))
  plot(gam_model, pages = 1, shade = TRUE, seWithMean = TRUE,
       main = pred_set_name)
  dev.off()
  
  # Model-metrics: AIC, BIC, Deviance, Deviance Explained,
  #    Splines (CWSI, SWC, VPD) Statistic + p-Wert, Tree.Species effects

  # Model-metrics
  mod_aic <- AIC(gam_model)
  mod_bic <- BIC(gam_model)
  
  # Deviance explained from summary
  mod_summary <- summary(gam_model)
  mod_dev_expl <- mod_summary$dev.expl  # Value also displayed in the summary
  
  # Spline-Statistics: Chi^2 and p-values
  sm <- mod_summary$s.table
  
  get_spline_stats <- function(sm_table, term_name) {
    if (!is.null(sm_table)) {
      idx <- grep(term_name, rownames(sm_table))
      if (length(idx) == 1) {
        chi_sq <- sm_table[idx, "Chi.sq"]
        pval   <- sm_table[idx, "p-value"]
      } else {
        chi_sq <- NA
        pval <- NA
      }
    } else {
      chi_sq <- NA
      pval <- NA
    }
    return(list(chi_sq = chi_sq, pval = pval))
  }
  
  s_CWSI <- get_spline_stats(sm, "CWSI")
  s_SWC  <- get_spline_stats(sm, "SWC")
  s_VPD  <- get_spline_stats(sm, "VPD")
  
  # Parametric coefficients for Tree.species
  coefs <- mod_summary$p.table
  
  get_coef_stats <- function(coef_table, term_name) {
    if (!is.null(coef_table)) {
      idx <- grep(term_name, rownames(coef_table))
      if (length(idx) == 1) {
        est <- coef_table[idx, "Estimate"]
        pval <- coef_table[idx, "Pr(>|z|)"]
      } else {
        est <- NA
        pval <- NA
      }
    } else {
      est <- NA
      pval <- NA
    }
    return(list(est = est, pval = pval))
  }
  
  tree_menz <- get_coef_stats(coefs, "Tree.SpeciesP. menziesii")
  tree_robur <- get_coef_stats(coefs, "Tree.SpeciesQ. robur")
  
  mod_r2adj <- as.numeric(sub(".*=\\s*|\\s+Deviance.*", "", 
                              grep("R-sq.*adj", capture.output(mod_summary), 
                                   value = TRUE)))
  
  
  # Save results in model_summary
  model_summary <- rbind(
    model_summary,
    data.frame(
      Predictor_Set = pred_set_name,
      AIC = mod_aic,
      BIC = mod_bic,
      R2adj = mod_summary$r.sq,
      Deviance_Explained = mod_dev_expl,
      s_CWSI_chi_sq = s_CWSI$chi_sq,
      s_CWSI_pvalue = s_CWSI$pval,
      s_SWC_chi_sq  = s_SWC$chi_sq,
      s_SWC_pvalue  = s_SWC$pval,
      s_VPD_chi_sq  = s_VPD$chi_sq,
      s_VPD_pvalue  = s_VPD$pval,
      Tree.SpeciesP.menziesii = tree_menz$est,
      Tree.SpeciesQ.robur = tree_robur$est,
      Tree.SpeciesP.menziesii_pval = tree_menz$pval,
      Tree.SpeciesQ.robur_pval = tree_robur$pval
    )
  )
  
  
  # Save model
  saveRDS(gam_model,
          file = paste0("./data/results/statistical_analysis/GAM_Model_",
                        pred_set_name, ".rds"))
  
}


# Save DHARMa Testresults
write.csv(dharma_summary, 
          "./data/results/statistical_analysis/DHARMa_Tests_GAM.csv", 
          row.names = FALSE)

write.csv(model_summary,
          "./data/results/statistical_analysis/DHARMa_GAM_model_summaries.csv",
          row.names = FALSE)

write.csv(oof_results_all_bal,
          "./data/results/statistical_analysis/DHARMa_GAM_oof_cv_results.csv",
          row.names = FALSE)


# write.csv(dharma_summary,
#           "./data/results/statistical_analysis_VPD_uncor/DHARMa_Tests_GAM.csv",
#           row.names = FALSE)
# 
# write.csv(model_summary,
#           "./data/results/statistical_analysis_VPD_uncor/DHARMa_GAM_model_summaries.csv",
#           row.names = FALSE)
# 
# write.csv(oof_results_all_bal,
#           "./data/results/statistical_analysis_VPD_uncor/DHARMa_GAM_oof_cv_results.csv",
#           row.names = FALSE)
# 
# write.csv(dharma_summary,
#           "./data/results/statistical_analysis_VPD_median_cor/DHARMa_Tests_GAM.csv",
#           row.names = FALSE)
# 
# write.csv(model_summary,
#           "./data/results/statistical_analysis_VPD_median_cor/DHARMa_GAM_model_summaries.csv",
#           row.names = FALSE)
# 
# write.csv(oof_results_all_bal,
#           "./data/results/statistical_analysis_VPD_median_cor/DHARMa_GAM_oof_cv_results.csv",
#           row.names = FALSE)


# Change predictor set names
predictor_names <- c(
  "CWSI_only" = "CWSI",
  "SWC_only" = "SWC",
  "VPD_only" = "VPD",
  "CWSI_SWC_VPD" = "CWSI+SWC+VPD",
  "CWSI_SWC" = "CWSI+SWC",
  "CWSI_VPD" = "CWSI+VPD",
  "SWC_VPD" = "SWC+VPD"
)

# format model summary
model_summary_formatted <- model_summary %>%
  mutate(Predictor_Set = as.character(Predictor_Set)) %>%   
  mutate(Predictor_Set = dplyr::recode(Predictor_Set, !!!predictor_names)) %>%
  arrange(desc(Deviance_Explained)) %>%
  mutate(
    AIC = round(AIC, 1),
    BIC = round(BIC, 1),
    R2adj = round(R2adj, 2),
    Deviance_Explained = round(Deviance_Explained, 2),
    s_CWSI_chi_sq = round(s_CWSI_chi_sq, 2),
    s_CWSI_pvalue = sapply(s_CWSI_pvalue, format_pval),
    s_SWC_chi_sq = round(s_SWC_chi_sq, 2),
    s_SWC_pvalue = sapply(s_SWC_pvalue, format_pval),
    s_VPD_chi_sq = round(s_VPD_chi_sq, 2),
    s_VPD_pvalue = sapply(s_VPD_pvalue, format_pval)
  )


# format oof_results_all_bal
oof_results_all_bal_GAM <- oof_results_all_bal %>%
  filter(Model == "GAM") %>%
  mutate(Predictor_Set = dplyr::recode(Predictor_Set, !!!predictor_names)) %>%
  dplyr::select(-c(Model, dist0N0, Tree.Species)) %>%
  dplyr::select(Balanced_Accuracy, everything()) %>%
  arrange(desc(Balanced_Accuracy)) %>%
  mutate(
    Balanced_Accuracy = round(Balanced_Accuracy, 4),
    Accuracy = round(Accuracy, 4),
    Kappa = round(Kappa, 4),
    Precision = round(Precision, 4),
    Recall = round(Recall, 4),
    F1 = round(F1, 4)
  )

# combine model_summary + oof_results 
GAM_results_combined <- left_join(oof_results_all_bal_GAM, 
                                  model_summary_formatted,
                                  by = "Predictor_Set")

# change to long format
GAM_long_df <- GAM_results_combined %>%
  mutate(across(-Predictor_Set, as.character)) %>%
  pivot_longer(cols = -Predictor_Set, names_to = "Metric", values_to = "Value")

# wide format
GAM_wide_df <- GAM_long_df %>%
  pivot_wider(names_from = Predictor_Set, values_from = Value)

# export csv
write.csv(GAM_results_combined, 
          "./data/results/statistical_analysis/GAM_results_combined.csv", 
          row.names = FALSE)
write.csv(GAM_wide_df, 
          "./data/results/statistical_analysis/GAM_results_wide.csv", 
          row.names = FALSE)


# write.csv(GAM_results_combined,
#           "./data/results/statistical_analysis_VPD_uncor/GAM_results_combined.csv",
#           row.names = FALSE)
# write.csv(GAM_wide_df,
#           "./data/results/statistical_analysis_VPD_uncor/GAM_results_wide.csv",
#           row.names = FALSE)
# 
# write.csv(GAM_results_combined,
#           "./data/results/statistical_analysis_VPD_median_cor/GAM_results_combined.csv",
#           row.names = FALSE)
# write.csv(GAM_wide_df,
#           "./data/results/statistical_analysis_VPD_median_cor/GAM_results_wide.csv",
#           row.names = FALSE)


# Calculate thresholds for P > 50% that TWDmin > 0 using GAM-CWSI

# create artifical values
cwsi_vals <- seq(0, 1, by = 0.001)
species_levels <- levels(data_B$Tree.Species)

gam_mod <- gam_list[[1]]
# empty list for all tree species
plots_list <- list()
thresholds_list <- list()

# loop over all tree species
for (sp in species_levels) {
  
  newdat <- data.frame(CWSI = cwsi_vals, Tree.Species = sp)
  
  # Prediction on Link-Scale
  pred <- predict(gam_mod, newdata = newdat, type = "link", se.fit = TRUE)
  
  # UChange to Response-Scale
  newdat <- newdat %>%
    mutate(
      fit_link = pred$fit,
      se_link  = pred$se.fit,
      Probability = plogis(fit_link),
      lower = plogis(fit_link - 1.96 * se_link),
      upper = plogis(fit_link + 1.96 * se_link)
    )
  
  # select threhsolds
  valid_idx <- which(!is.na(newdat$Probability))
  prob_valid <- newdat$Probability[valid_idx]
  cwsi_valid <- newdat$CWSI[valid_idx]
  
  if (any(prob_valid >= 0.5) && any(prob_valid <= 0.5)) {
    cwsi_thr_mid <- approx(x = prob_valid, y = cwsi_valid, xout = 0.5)$y
    ci_idx <- which(newdat$upper >= 0.5 & newdat$lower <= 0.5)
    cwsi_thr_min <- min(newdat$CWSI[ci_idx])
    cwsi_thr_max <- max(newdat$CWSI[ci_idx])
    
    newdat <- newdat %>%
      mutate(
        lower_thr = ifelse(CWSI >= cwsi_thr_min & CWSI <= cwsi_thr_max, lower, NA),
        upper_thr = ifelse(CWSI >= cwsi_thr_min & CWSI <= cwsi_thr_max, upper, NA)
      )
    
  } else {
    cwsi_thr_mid <- NA
    cwsi_thr_min <- NA
    cwsi_thr_max <- NA
    newdat$lower_thr <- NA
    newdat$upper_thr <- NA
  }
  
  thresholds_list[[sp]] <- data.frame(
    Tree.Species = sp,
    thr_mid = cwsi_thr_mid,
    thr_min = cwsi_thr_min,
    thr_max = cwsi_thr_max
  )
  
  # Create threshold plot per tree species
  p_sp <- ggplot(newdat, aes(x = CWSI, y = Probability)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "red", alpha = 0.2) +
    geom_line(color = "red", size = 1.5) +
    annotate("rect",
             xmin = cwsi_thr_min,
             xmax = cwsi_thr_max,
             ymin = -Inf,
             ymax = Inf,
             fill = "lightblue",
             alpha = 0.4) +
    geom_vline(xintercept = cwsi_thr_mid, linetype = "dashed", 
               color = "blue", linewidth = 1.2) +
    geom_hline(yintercept = 0.5, linetype = "dashed", 
               color = "grey40", linewidth = 1) +
    labs(
      title = paste0("Predicted probability for ", sp),
      x = "CWSI",
      y = "Probability (TWDmin > 0)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      axis.text = element_text(size = 13)
    )
  
  plots_list[[sp]] <- p_sp
}

# Collect thresholds in overview table
thresholds_df <- bind_rows(thresholds_list)
print(thresholds_df)

# Combined Plot using all tree species (facetting)
# Collect data from all tree species
all_newdat <- do.call(rbind, lapply(names(plots_list), function(sp) {
  df <- plots_list[[sp]]$data
  df$Tree.Species <- sp
  df
}))

species_colors <- c(
  "F. sylvatica" = "#298c8c",
  "P. menziesii" = "orange",
  "Q. robur" = "darkmagenta"
)

# Order thresholds labeling
thresholds_df <- thresholds_df %>%
  arrange(desc(thr_mid)) %>% 
  mutate(
    y_pos = seq(0.15, 0.35, length.out = n())  # define distance between texts
  )

# Threshold plot (Figure 5)
gg_combined <- ggplot(all_newdat, 
                      aes(x = CWSI, y = Probability, color = Tree.Species)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Tree.Species), 
              alpha = 0.2, color = NA) +
  geom_hline(yintercept = 0.5, color = "grey40", size = 1) +
  geom_vline(data = thresholds_df, 
             aes(xintercept = thr_mid, color = Tree.Species),
             linetype = "dashed", size = 1.2, show.legend = FALSE ) +
  
  geom_text(
    data = thresholds_df,
    aes(x = 0.7, y = y_pos, 
        label = paste0("Thr = ", sprintf("%.2f", thr_mid)), 
        color = Tree.Species),
    angle = 0,
    size = 7,
    hjust = 0,          
    inherit.aes = FALSE, 
    show.legend = FALSE,
    fontface = "bold" 
  ) +
  
  labs(
    x = "CWSI",
    y = expression(bold("P ("*TWD[min]*" > 0)")),
    color = "Tree species",
    fill = "Tree species"
  ) +
  
  scale_color_manual(values = species_colors) +
  scale_fill_manual(values = species_colors) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 13)
  )

print(gg_combined)

ggsave("./graphics/results/thresholds/CWSI_irrig_thres_GAM_all_species.png",
       plot = gg_combined, width = 12, height = 7, dpi = 300)



# Calculate thresholds for P > 50% that TWDmin > 0 using GAM-SWC

# Create artifical values
swc_vals <- seq(0,50, length.out = 1000) 
species_levels <- levels(data_B$Tree.Species)

gam_mod <- gam_list[[5]]
plots_list_swc <- list()
thresholds_list_swc <- list()

# Loop over all tree species
for (sp in species_levels) {
  
  newdat <- data.frame(SWC = swc_vals, Tree.Species = sp)
  
  # Predictions on Link-Scale
  pred <- predict(gam_mod, newdata = newdat, type = "link", se.fit = TRUE)
  
  # Change to Response-Scale
  newdat <- newdat %>%
    mutate(
      fit_link = pred$fit,
      se_link  = pred$se.fit,
      Probability = plogis(fit_link),
      lower = plogis(fit_link - 1.96 * se_link),
      upper = plogis(fit_link + 1.96 * se_link)
    )
  
  # Define thresholds
  valid_idx <- which(!is.na(newdat$Probability))
  prob_valid <- newdat$Probability[valid_idx]
  swc_valid <- newdat$SWC[valid_idx]
  
  if (any(prob_valid >= 0.5) && any(prob_valid <= 0.5)) {
    swc_thr_mid <- approx(x = prob_valid, y = swc_valid, xout = 0.5)$y
    ci_idx <- which(newdat$upper >= 0.5 & newdat$lower <= 0.5)
    swc_thr_min <- min(newdat$SWC[ci_idx])
    swc_thr_max <- max(newdat$SWC[ci_idx])
    
    newdat <- newdat %>%
      mutate(
        lower_thr = ifelse(SWC >= swc_thr_min & SWC <= swc_thr_max, lower, NA),
        upper_thr = ifelse(SWC >= swc_thr_min & SWC <= swc_thr_max, upper, NA)
      )
    
  } else {
    swc_thr_mid <- NA
    swc_thr_min <- NA
    swc_thr_max <- NA
    newdat$lower_thr <- NA
    newdat$upper_thr <- NA
  }
  
  thresholds_list_swc[[sp]] <- data.frame(
    Tree.Species = sp,
    thr_mid = swc_thr_mid,
    thr_min = swc_thr_min,
    thr_max = swc_thr_max
  )
  
  # Plot per tree species with threshold
  p_sp <- ggplot(newdat, aes(x = SWC, y = Probability)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "red", alpha = 0.2) +
    geom_line(color = "red", size = 1.5) +
    annotate("rect",
             xmin = swc_thr_min,
             xmax = swc_thr_max,
             ymin = -Inf,
             ymax = Inf,
             fill = "lightblue",
             alpha = 0.4) +
    geom_vline(xintercept = swc_thr_mid, linetype = "dashed", 
               color = "blue", size = 1.2) +
    geom_hline(yintercept = 0.5, linetype = "dashed", 
               color = "grey40", size = 1) +
    labs(
      title = paste0("Predicted probability for ", sp),
      x = "SWC",
      y = "Probability (TWDmin > 0)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      axis.text = element_text(size = 13)
    )
  
  plots_list_swc[[sp]] <- p_sp
}

# Create overview table with thresholds
thresholds_swc_df <- bind_rows(thresholds_list_swc)
print(thresholds_swc_df)

# Combined plot using all tree species
all_newdat_swc <- do.call(rbind, lapply(names(plots_list_swc), function(sp) {
  df <- plots_list_swc[[sp]]$data
  df$Tree.Species <- sp
  df
}))

# Arrange threshold-labeling
thresholds_swc_df <- thresholds_swc_df %>%
  arrange(desc(thr_mid)) %>%
  mutate(y_pos = seq(0.15, 0.35, length.out = n()))

# Threshold plot (S 24)
gg_combined_swc <- ggplot(all_newdat_swc, 
                          aes(x = SWC, y = Probability, color = Tree.Species)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Tree.Species), 
              alpha = 0.2, color = NA) +
  geom_hline(yintercept = 0.5, color = "grey40", size = 1) +
  geom_vline(data = thresholds_swc_df, aes(xintercept = thr_mid, 
                                           color = Tree.Species),
             linetype = "dashed", size = 1.2, show.legend = FALSE) +
  geom_text(
    data = thresholds_swc_df,
    aes(x = 26, y = y_pos + 0.5,
        label = paste0("Thr = ", sprintf("%.2f", thr_mid)), 
        color = Tree.Species),
    angle = 0, size = 7, hjust = 0, inherit.aes = FALSE,
    show.legend = FALSE, fontface = "bold"
  ) +
  labs(
    x = "SWC [%]",
    y = expression(bold("P ("*TWD[min]*" > 0)")),
    color = "Tree species",
    fill = "Tree species"
  ) +
  scale_color_manual(values = species_colors) +
  scale_fill_manual(values = species_colors) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 13)
  )

print(gg_combined_swc)

ggsave("./graphics/results/thresholds/SWC_irrig_thres_GAM_all_species.png",
       plot = gg_combined_swc, width = 12, height = 7, dpi = 300)


# Calculate thresholds for P > 50% that TWDmin > 0 using GAM-CWSI+SWC+VPD

# Create artificial prediction data
contour_df <- expand.grid(
  CWSI = c(0.1, 0.5, 0.9),
  SWC  = seq(0, 50, length.out = 200),
  VPD_kPa = seq(0, 4, length.out = 40),
  Tree.Species = levels(data_B$Tree.Species)
)

# Predict
contour_df$Probability <- predict(gam_list[[4]], newdata = contour_df, 
                                  type = "response")

# Line at 0.5 ± tolerance
contour_lines <- contour_df %>%
  group_by(CWSI, Tree.Species) %>%
  filter(abs(Probability - 0.5) < 0.01) %>%
  ungroup() %>%
  mutate(
    CWSI_label = case_when(
      CWSI ==  0.1 ~ "0.1",
      CWSI == 0.5  ~ "0.5",
      CWSI == 0.9 ~ "0.9"
    ),
    linetype_value = case_when(
      CWSI_label == "0.5" ~ "solid",
      CWSI_label == "0.1" ~ "dashed",
      CWSI_label == "0.9" ~ "dotted"
    )
  )

# Create plot (Figure 6) with different linetypes and sizes per CWSI value
ggplot(contour_lines, aes(x = VPD_kPa, y = SWC, color = Tree.Species)) +
  geom_line(aes(linetype = CWSI_label, size = CWSI_label)) +
  scale_color_manual(values = species_colors) +
  scale_linetype_manual(values = c("0.1" = "dashed", 
                                   "0.5" = "solid", 
                                   "0.9" = "dotdash")) +
  scale_size_manual(values = c("0.1" = 0.6, 
                               "0.5" = 1.2, 
                               "0.9" = 0.6)) +
  labs(
    x = "VPD [kPa]",
    y = "SWC [%]",
    color = "Tree Species",
    linetype = "CWSI",
    size = "CWSI"
  ) +
  xlim(0.5, 4)+
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5)
  )

ggsave("./graphics/results/thresholds/ThrLines_SWCVPD_CWSI_values_species_Linetype.png",
       width = 10, height = 6, dpi = 300)

# ggsave("./graphics/results/thresholds/ThrLines_SWCVPD_CWSI_values_species_Linetype_VPD_uncor.png",
#        width = 10, height = 6, dpi = 300)
# 
# ggsave("./graphics/results/thresholds/ThrLines_SWCVPD_CWSI_values_species_Linetype_VPD_median_cor.png",
#        width = 10, height = 6, dpi = 300)

# Calculate thresholds for P > 50% that TWDmin > 0 using GAM-CWSI+SWC
mod <- gam_list[[2]]

# Create grid for heatmap 
grid_df <- expand.grid(
  CWSI = seq(0, 1, length.out = 100),
  SWC  = seq(0, 50, length.out = 500),
  Tree.Species = levels(data_B$Tree.Species)
)

# Predict
grid_df$Probability <- predict(gam_list[[2]], newdata = grid_df, 
                               type = "response")


# Plot S 26: P= 0.5-Thresholds + confidence band
ggplot(grid_df, aes(x = CWSI, y = SWC)) +
  # confidence band as shaded area between P=0.475 and Prob=0.525 
  # -> 5 % Band
  geom_contour_filled(
    aes(z = Probability, fill = Tree.Species),
    breaks = c(0.5 - 0.025, 0.5 + 0.025), 
    alpha = 0.2
  ) +
  # 0.5-Line
  geom_contour(
    aes(z = Probability, color = Tree.Species),
    breaks = 0.5,
    size = 1
  ) +
  scale_color_manual(values = species_colors) +
  scale_fill_manual(values = species_colors) +
  labs(
    x = "CWSI",
    y = "SWC [%]",
    color = "Tree Species",
    fill = "Tree Species"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 14, face = "bold"),)

ggsave("./graphics/results/thresholds/ThrLines+perc_bands_SWCCWSI.png",
       width = 10, height = 6, dpi = 300)



# calculate cwsi means and sd betweeen TWDmin classes
CWSI_summary_table_classes <- data_B %>%
  mutate(class = dplyr::recode(class, "nonzero" = "> 0", "zero" = "0")) %>%
  group_by(Tree.Species, class) %>%
  summarise(
    mean_CWSI = mean(CWSI, na.rm = TRUE),
    sd_CWSI = sd(CWSI, na.rm = TRUE),
    n = n()
  ) %>%
  ungroup()

CWSI_summary_table_classes
write.csv(CWSI_summary_table_classes, 
          "./data/results/statistical_analysis/CWSI_summary_table_classes.csv")

CWSI_summary_table_classes <- data_B %>%
  mutate(class = dplyr::recode(class, "nonzero" = "> 0", "zero" = "0")) %>%
  group_by(Tree.Species, class) %>%
  summarise(
    mean_CWSI = mean(CWSI, na.rm = TRUE),
    sd_CWSI = sd(CWSI, na.rm = TRUE),
    n = n()
  ) %>%
  ungroup()

# Calculate differences: "TWDmin > 0" minus "TWDmin = 0"
CWSI_summary_diff <- CWSI_summary_table_classes %>%
  pivot_wider(
    names_from = class,
    values_from = c(mean_CWSI, sd_CWSI, n),
    names_sep = "_"
  ) %>%
  mutate(
    diff_mean = `mean_CWSI_> 0` - `mean_CWSI_0`,
    mean_sd_0 = paste0(round(`mean_CWSI_0`, 3), " ± ", round(`sd_CWSI_0`, 3)),
    mean_sd_nonzero = paste0(round(`mean_CWSI_> 0`, 3), " ± ", round(`sd_CWSI_> 0`, 3))
  )

# recalculate differences values in %
CWSI_summary_diff <- CWSI_summary_diff %>%
  mutate(diff_mean_perc = diff_mean/mean_CWSI_0)

write.csv(CWSI_summary_diff, 
          "./data/results/statistical_analysis/CWSI_summary_table_classes_difs.csv")


# Plot differences in CWSI, SWC and VPD between MinTWD = 0 and MinTWD > 0
# define colors
cols <- c("   0" = "blue", "> 0" = "red")

# Boxplot with Wilcoxon-Test ( **** p < 0.0001, ** < 0.01) (Figure 4) to
# show significant CWSI differences between classes
p <- ggplot(data_B %>%
              mutate(class = dplyr::recode(class, "nonzero" = "> 0"),
                     class = dplyr::recode(class, "zero" = "   0")), 
            aes(x = Tree.Species, y = CWSI, fill = class)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.size = 2, 
               color = "black") +
  
  # Wilcoxon-test per tree species
  stat_compare_means(
    aes(group = class),
    method = "wilcox.test",
    label = "p.signif",        # shows “ns”, “*”, “**”, “***”
    label.y = max(data_B$CWSI, na.rm = TRUE) * 1.05, # position above boxes
    hide.ns = TRUE
  ) +

  labs(
    title = "CWSI per tree species and TWDmin-class",
    x = "Tree species",
    y = "CWSI",
    fill = expression(TWD[min])  
  ) +
  scale_fill_manual(values = cols) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 13), # horizontal
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  )

print(p)

ggsave("./graphics/results/CWSI_difs_MinTWD0N0.png",
       plot = p, width = 10, height = 5, dpi = 300)

# Boxplot with Wilcoxon-Test ( **** p < 0.0001, ** < 0.01) for SWC differences
p <- ggplot(data_B %>%
              mutate(class = dplyr::recode(class, "nonzero" = "> 0"),
                     class = dplyr::recode(class, "zero" = "   0")), 
            aes(x = Tree.Species, y = SWC, fill = class)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.size = 2, 
               color = "black") +
  
  # Wilcoxon-test per tree species
  stat_compare_means(
    aes(group = class),
    method = "wilcox.test",
    label = "p.signif",       
    label.y = max(data_B$SWC, na.rm = TRUE) * 1.05, 
    hide.ns = TRUE
  ) +

  labs(
    title = "SWC per tree species and TWDmin-class",
    x = "Tree species",
    y = "SWC [%]",
    fill = expression(TWD[min])  
  ) +
  scale_fill_manual(values = cols) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 13), # horizontal
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  )

print(p)

ggsave("./graphics/results/SWC_difs_MinTWD0N0.png",
       plot = p, width = 10, height = 5, dpi = 300)

# Boxplot with Wilcoxon-Test ( **** p < 0.0001, ** < 0.01) for VPD differences
p <- ggplot(data_B %>%
              mutate(class = dplyr::recode(class, "nonzero" = "> 0"),
                     class = dplyr::recode(class, "zero" = "   0")), 
            aes(x = Tree.Species, y = VPD_kPa, fill = class)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.size = 2, 
               color = "black") +
  
  # Wilcoxon-Test per tree species
  stat_compare_means(
    aes(group = class),
    method = "wilcox.test",
    label = "p.signif",        
    label.y = max(data_B$VPD_kPa, na.rm = TRUE) * 1.05, 
    hide.ns = TRUE
  ) +
  labs(
    title = "VPD per tree species and TWDmin-class",
    x = "Tree species",
    y = "VPD [kPa]",
    fill = expression(TWD[min])  
  ) +
  scale_fill_manual(values = cols) +

  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 13), 
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  )

print(p)

ggsave("./graphics/results/VPD_difs_MinTWD0N0.png",
       plot = p, width = 10, height = 5, dpi = 300)

# ggsave("./graphics/results/VPD_difs_MinTWD0N0_VPD_uncor.png",
#        plot = p, width = 10, height = 5, dpi = 300)
# 
# ggsave("./graphics/results/VPD_difs_MinTWD0N0_VPD_median_cor.png",
#        plot = p, width = 10, height = 5, dpi = 300)

# Combined plot with TWDmin class differences for CWSI, SWC and VPD
data_B_mod <- data_B %>%
  mutate(class = dplyr::recode(class,
                        "nonzero" = "> 0",
                        "zero" = "   0"))


# Create single plots
p1 <- make_boxplot(data_B_mod, "CWSI", "CWSI", show_x = FALSE)
p2 <- make_boxplot(data_B_mod, "SWC", "SWC [%]", show_x = FALSE)
p3 <- make_boxplot(data_B_mod, "VPD_kPa", "VPD [kPa]", show_x = TRUE)  

# Combine plots with minimal spacer inbetween
combined_plot <- p1 / plot_spacer() / p2 / plot_spacer() / p3 +
  plot_layout(guides = "collect", 
              heights = c(1, 0.05, 1, 0.05, 1)) &  # small spacer 
  theme(legend.position = "right",
        legend.box = "vertical")

print(combined_plot)

ggsave("./graphics/results/CWSI_SWC_VPD_difs_MinTWD0N0.png",
       plot = combined_plot, width = 12, height = 8, dpi = 300)

# Figure S 24: SWC and VPD differences according to TWDmin classes
combined_plot2 <- p2 / plot_spacer() / p3 +
  plot_layout(guides = "collect", heights = c(1, 0.05, 1, 0.05, 1)) &  
  theme(legend.position = "right",
        legend.box = "vertical")

print(combined_plot2)

ggsave("./graphics/results/SWC_VPD_difs_MinTWD0N0.png",
       plot = combined_plot2, width = 12, height = 8, dpi = 300)

