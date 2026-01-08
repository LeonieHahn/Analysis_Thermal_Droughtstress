# - k-means clustering to differentiate reference panels from background if 
#   there is still background (e. g. white) included
# - remove values of one of the mean-clusters only if the differences between 
#   the means of the kmean clusters is larger than 0.5
# - use colder cluster for green_wet and warmer cluster for green_dry
# - remove values outside 10./90. percentile range (=outliers/background)
# - compare the usage of other quantiles
# - create dataframe with statistics of temperature values
# - check the preprocessed images 
#   (after k-means clustering & percentile clipping) and compare to originals

library(terra)
library(fs)
library(dplyr)
library(readr)
library(ggplot2)
library(viridisLite)
library(stringr)

source("./Code/Functions.R")

input_root <- "D:/Bewaesserung_Forstkulturen/Daten/Gewaechshaus/Gewaechshaus_Trockenstress2023/Thermal/GH_Analyse_Thermal/Cropped_Thermals_1eroded_pixel"

output_root <- "D:/Bewaesserung_Forstkulturen/Daten/Gewaechshaus/Gewaechshaus_Trockenstress2023/Thermal/GH_Analyse_Thermal/Cropped_Thermals_1eroded_pixel_Refs_kmeans"

cluster_means_df <- data.frame(
  file = character(),
  date = character(),
  mean_cluster1 = numeric(),
  mean_cluster2 = numeric(),
  stringsAsFactors = FALSE
)

subdirs <- dir_ls(input_root, type = "directory", recurse = FALSE)

# Loop through all folders and perform k-means clustering for the reference 
# surfaces 
for (subdir in subdirs) {
  # extract date from folder
  date <- path_file(subdir)
  
  # create corresponding folder for results
  out_subdir <- path(output_root, date)
  dir_create(out_subdir)
  
  # list all TIF-files in the subfolder
  tif_files <- dir_ls(subdir, glob = "*.tif")
  
  for (tif_path in tif_files) {
    message("Process: ", tif_path)
    
    fname <- path_file(tif_path)
    
    # Check, if the filename corresponds to the object type 
    # that should be processed
    if (!grepl("Green_Wet_Shade|Green_Wet_Sun|Green_Dry_Shade|Green_Dry_Sun", 
               fname)) {
      message("No relevant object, original image is copied")
      out_file <- path(out_subdir, fname)
      file_copy(tif_path, out_file, overwrite = TRUE)
      next
    }
    
    img <- rast(tif_path)
    vals <- values(img)
    valid_idx <- which(!is.na(vals))
    valid_vals <- vals[valid_idx, , drop = FALSE]
    
    if (nrow(valid_vals) < 2) {
      message("Not enough valid pixels, skip...")
      next
    }
    
    set.seed(42)
    km <- kmeans(valid_vals, centers = 2)
    
    mean_cl1 <- mean(valid_vals[km$cluster == 1])
    mean_cl2 <- mean(valid_vals[km$cluster == 2])
    
    cluster_means_df <- rbind(cluster_means_df, data.frame(
      file = fname,
      date = date,
      mean_cluster1 = mean_cl1,
      mean_cluster2 = mean_cl2,
      stringsAsFactors = FALSE
    ))
    
    # check, if differences between means of cluster are > 0.5
    if (abs(mean_cl1 - mean_cl2) > 0.5) {
      # keep colder cluster
      if (grepl("Green_Wet_Shade|Green_Wet_Sun", fname)) {
        tree_cluster <- if (mean_cl1 < mean_cl2) 1 else 2
      } else {
        # keep warmer cluster for green_dry
        tree_cluster <- if (mean_cl1 > mean_cl2) 1 else 2
      }
      
      mask_vals <- rep(NA, length(vals))
      mask_vals[valid_idx[km$cluster == tree_cluster]] <- 
        vals[valid_idx[km$cluster == tree_cluster]]
      
      tree_img <- img
      values(tree_img) <- mask_vals
      out_file <- path(out_subdir, fname)
      writeRaster(tree_img, out_file, overwrite = TRUE)
      
    } else {
      out_file <- path(out_subdir, fname)
      writeRaster(img, out_file, overwrite = TRUE)
    }
  }}


cluster_means_df <- cluster_means_df %>%
  mutate(mean_dif = mean_cluster1-mean_cluster2,
         mean_dif_abs = abs(mean_dif))

cluster_means_df <- cluster_means_df %>%
  mutate(
    Type = str_extract(file, "Tree|Green_Wet|Green_Dry"),
    Position = str_extract(file, "Shade|Sun")
  )

cluster_means_df_Tree <- cluster_means_df %>%
  filter(Type == "Tree")

cluster_means_stats <- cluster_means_df %>%
  group_by(Type) %>%
  summarize(mean_difs = mean(mean_dif_abs),
            median_difs = median(mean_dif_abs))


# Input and output folders
input_root <- "D:/Bewaesserung_Forstkulturen/Daten/Gewaechshaus/Gewaechshaus_Trockenstress2023/Thermal/GH_Analyse_Thermal/Cropped_Thermals_1eroded_pixel_Refs_kmeans"
output_root <- "D:/Bewaesserung_Forstkulturen/Daten/Gewaechshaus/Gewaechshaus_Trockenstress2023/Thermal/GH_Analyse_Thermal/Cropped_Thermals_percentile_clipped"

# Define all quantile ranges that should be tested
quantile_ranges <- list(
  Q10_Q90   = c(0.10, 0.90),
  Q5_Q95    = c(0.05, 0.95),
  Q1_Q99    = c(0.01, 0.99),
  Q15_Q85   = c(0.15, 0.85),
  IQR_based = "IQR"
)

# Empty result dataframe
percentile_means_df <- data.frame(
  file = character(),
  Date = character(),
  object = character(),
  Position = character(),
  Tree.ID = character(),
  Treatment = character(),
  Tree.Species = character(),
  RangeType = character(),
  mean_objectarea = numeric(),
  median_objectarea = numeric(),
  min_objectarea = numeric(),
  max_objectarea = numeric(),
  sd_objectarea = numeric(),
  pixelamount_used = numeric(),
  stringsAsFactors = FALSE
)

# Process all subdirectories
subdirs <- dir_ls(input_root, type = "directory", recurse = FALSE)

for (subdir in subdirs) {
  folder_name <- path_file(subdir)
  Date <- substr(folder_name, 1, 8)
  
  object_part <- sub("^[0-9]+_", "", folder_name)
  object <- sub("_(Shade|Sun)$", "", object_part)
  position <- ifelse(grepl("Shade", folder_name), "Shade", "Sun")
  
  tif_files <- dir_ls(subdir, glob = "*.tif")
  
  for (tif_path in tif_files) {
    message("Process: ", tif_path)
    img <- rast(tif_path)
    vals <- values(img)
    valid_idx <- which(!is.na(vals))
    valid_vals <- vals[valid_idx]
    
    if (length(valid_vals) < 10) {
      message("Not enough final pixels, skip...")
      next
    }
    
    # Tree-infos from filename
    filename <- path_file(tif_path)
    id_match <- regmatches(filename, 
                           regexec("crop_([A-Z][0-9]+)_([DW])", filename))[[1]]
    tree_id <- paste0(id_match[2], "_", id_match[3])
    behandlung <- id_match[3]
    art_code <- substr(id_match[2], 1, 1)
    baumart <- dplyr::case_when(
      art_code == "B" ~ "Fagus sylvatica",
      art_code == "T" ~ "Abies alba",
      art_code == "E" ~ "Quercus robur",
      art_code == "D" ~ "Pseudotsuga menziesii",
      TRUE ~ NA_character_
    )
    
    # Process every range
    for (range_name in names(quantile_ranges)) {
      mask_vals <- rep(NA, length(vals))
      
      if (range_name == "IQR_based") {
        q25 <- quantile(valid_vals, 0.25, na.rm = TRUE)
        q75 <- quantile(valid_vals, 0.75, na.rm = TRUE)
        iqr <- q75 - q25
        lower <- q25 - 1.5 * iqr
        upper <- q75 + 1.5 * iqr
      } else {
        lower <- quantile(valid_vals, quantile_ranges[[range_name]][1], 
                          na.rm = TRUE)
        upper <- quantile(valid_vals, quantile_ranges[[range_name]][2], 
                          na.rm = TRUE)
      }
      
      # Clipping-rules:
      # - remove upper (warmer) values for the wet reference surfaces
      # - remove lower (colder) values for the dry reference surfaces
      # - remove lower (colder) and upper (warmer) values for the tree saplings
      if (object_part %in% c("Green_Wet_Shade", "Green_Wet_Sun")) {
        keep_idx <- valid_idx[valid_vals <= upper]
      } else if (object_part %in% c("Green_Dry_Shade", "Green_Dry_Sun")) {
        keep_idx <- valid_idx[valid_vals >= lower]
      } else if (object_part %in% c("Tree_Shade", "Tree_Sun")) {
        keep_idx <- valid_idx[valid_vals >= lower & valid_vals <= upper]
      } else {
        keep_idx <- valid_idx
      }
      
      mask_vals[keep_idx] <- vals[keep_idx]
      
      # Calculate statistical metrics for the pixels extracted per object
      mean_object <- mean(mask_vals, na.rm = TRUE)
      median_object <- median(mask_vals, na.rm = TRUE)
      min_object <- min(mask_vals, na.rm = TRUE)
      max_object <- max(mask_vals, na.rm = TRUE)
      sd_object <- sd(mask_vals, na.rm = TRUE)
      pix_amount <- length(keep_idx)
      
      # Save results
      percentile_means_df <- rbind(percentile_means_df, data.frame(
        file = filename,
        Date = Date,
        object = object,
        Position = position,
        Tree.ID = tree_id,
        Treatment = behandlung,
        Tree.Species = baumart,
        RangeType = range_name,
        mean_objectarea = mean_object,
        median_objectarea = median_object,
        min_objectarea = min_object,
        max_objectarea = max_object,
        sd_objectarea = sd_object,
        pixelamount_used = pix_amount,
        stringsAsFactors = FALSE
      ))
      
      # Save image in folder according to percentile range
      out_subdir <- path(output_root, paste0(folder_name, "_", range_name))
      dir_create(out_subdir)
      clipped_img <- img
      values(clipped_img) <- mask_vals
      out_file <- path(out_subdir, filename)
      writeRaster(clipped_img, out_file, overwrite = TRUE)
    }
  }
}

message("Processing done.")

# Postprocessing of dataframe with results
percentile_means_df <- percentile_means_df %>%
  mutate(mean_median_dif = mean_objectarea - median_objectarea,
         range = max_objectarea - min_objectarea,
         Date = as.Date(as.character(Date), format = "%Y%m%d"))

write.csv(percentile_means_df, 
          file = "./data/Thermal/results/percentile_multi_clipped_stats.csv", 
          row.names = FALSE)

percentile_means_df$Date <- as.Date(percentile_means_df$Date, format = "%Y%m%d")

# count how many pixels were used after removing the values outside the 
# 10% and 90% quantile
pixel_amount_used <- percentile_means_df %>%
  filter(RangeType == "Q10_Q90") %>%
  group_by(object) %>%
  summarise(
    mean_pixels = mean(pixelamount_used, na.rm = TRUE),
    min_pixels  = min(pixelamount_used, na.rm = TRUE),
    max_pixels  = max(pixelamount_used, na.rm = TRUE),
    .groups = "drop"
  )
  

# Add median values of the thermal images including the masked versions without 
# postprocessing

allpixels_root <- "D:/Bewaesserung_Forstkulturen/Daten/Gewaechshaus/Gewaechshaus_Trockenstress2023/Thermal/GH_Analyse_Thermal/Cropped_Thermals_all_pixels"

# use all subfolders as earlier
subdirs_allpixels <- dir_ls(allpixels_root, type = "directory", recurse = FALSE)

# Create empty dataframe for median values
allpixels_df <- data.frame(
  file = character(),
  Date = character(),
  object = character(),
  Position = character(),
  Tree.ID = character(),
  Treatment = character(),
  Tree.Species = character(),
  RangeType = character(),
  mean_objectarea = numeric(),
  median_objectarea = numeric(),
  min_objectarea = numeric(),
  max_objectarea = numeric(),
  sd_objectarea = numeric(),
  pixelamount_used = numeric(),
  stringsAsFactors = FALSE
)

for (subdir in subdirs_allpixels) {
  folder_name <- path_file(subdir)
  Date <- substr(folder_name, 1, 8)
  
  object_part <- sub("^[0-9]+_", "", folder_name)
  object <- sub("_(Shade|Sun)$", "", object_part)
  position <- ifelse(grepl("Shade", folder_name), "Shade", "Sun")
  
  tif_files <- dir_ls(subdir, glob = "*.tif")
  
  for (tif_path in tif_files) {
    message("Calculate median for: ", tif_path)
    img <- rast(tif_path)
    vals <- values(img)
    valid_vals <- vals[!is.na(vals)]
    
    if (length(valid_vals) < 10) next
    
    # Extract tree info from filename
    filename <- path_file(tif_path)
    id_match <- regmatches(filename, 
                           regexec("crop_([A-Z][0-9]+)_([DW])", filename))[[1]]
    tree_id <- paste0(id_match[2], "_", id_match[3])
    behandlung <- id_match[3]
    art_code <- substr(id_match[2], 1, 1)
    baumart <- dplyr::case_when(
      art_code == "B" ~ "Fagus sylvatica",
      art_code == "T" ~ "Abies alba",
      art_code == "E" ~ "Quercus robur",
      art_code == "D" ~ "Pseudotsuga menziesii",
      TRUE ~ NA_character_
    )
    
    # Calculate statistical metrics for values without percentile clipping
    mean_object <- mean(valid_vals, na.rm = TRUE)
    median_object <- median(valid_vals, na.rm = TRUE)
    min_object <- min(valid_vals, na.rm = TRUE)
    max_object <- max(valid_vals, na.rm = TRUE)
    sd_object <- sd(valid_vals, na.rm = TRUE)
    pix_amount <- length(valid_vals)
    
    # Save results
    allpixels_df <- rbind(allpixels_df, data.frame(
      file = filename,
      Date = Date,
      object = object,
      Position = position,
      Tree.ID = tree_id,
      Treatment = behandlung,
      Tree.Species = baumart,
      RangeType = "All_Pixels",
      mean_objectarea = mean_object,
      median_objectarea = median_object,
      min_objectarea = min_object,
      max_objectarea = max_object,
      sd_objectarea = sd_object,
      pixelamount_used = pix_amount,
      stringsAsFactors = FALSE
    ))
  }
}

# format date
allpixels_df$Date <- as.Date(as.character(allpixels_df$Date), format = "%Y%m%d")

# merge with existing table

# load results
combined_df <- bind_rows(percentile_means_df, allpixels_df)

# Calculate mean-median difference and value range
combined_df <- combined_df %>%
  mutate(mean_median_dif = mean_objectarea - median_objectarea,
         range = max_objectarea - min_objectarea)

# Save total result table
write.csv(combined_df,
          file = "./data/Thermal/results/percentile_and_allpixels_stats.csv",
          row.names = FALSE)


# view results and compare with original cropped thermal image
input_root <- "D:/Bewaesserung_Forstkulturen/Daten/Gewaechshaus/Gewaechshaus_Trockenstress2023/Thermal/GH_Analyse_Thermal/Cropped_Thermals_all_pixels"
output_root <- "D:/Bewaesserung_Forstkulturen/Daten/Gewaechshaus/Gewaechshaus_Trockenstress2023/Thermal/GH_Analyse_Thermal/Cropped_Thermals_percentile_clipped"

target_root <- "D:/Bewaesserung_Forstkulturen/Daten/Gewaechshaus/Gewaechshaus_Trockenstress2023/Thermal/GH_Analyse_Thermal/Cropped_Thermals_1eroded_pixel_Refs_kmeans_percentile1090_clipped"

# Copy all Q10_Q90 folders in new folder for comparison
all_subfolders <- list.dirs(output_root, recursive = FALSE, full.names = TRUE)
q10_folders <- all_subfolders[grepl("Q10_Q90", all_subfolders)]

# Create new directory
if(!dir.exists(target_root)) dir.create(target_root, recursive = TRUE)

# Copy all Q10_Q90-folder and rename
for(folder in q10_folders) {
  # Rename folder directory, so "_Q10_Q90" is removed
  new_name <- gsub("_Q10_Q90", "", basename(folder))
  dest_folder <- file.path(target_root, new_name)
  
  cat("Copying", folder, "to", dest_folder, "\n")
  dir.create(dest_folder, recursive = TRUE, showWarnings = FALSE)
  
  #  Copy all files in folder
  files_to_copy <- list.files(folder, full.names = TRUE)
  file.copy(files_to_copy, dest_folder, overwrite = TRUE)
}

cat("Finished copying all Q10_Q90 folders.\n")


input_files <- dir(input_root, pattern = "\\.tif$", recursive = TRUE, 
                   full.names = TRUE)

input_root_kmeans <- "D:/Bewaesserung_Forstkulturen/Daten/Gewaechshaus/Gewaechshaus_Trockenstress2023/Thermal/GH_Analyse_Thermal/Cropped_Thermals_1eroded_pixel_Refs_kmeans"

input_files_kmeans <- dir(input_root, pattern = "\\.tif$", recursive = TRUE, full.names = TRUE)



# compare k-means clustered image with its original version
start_index <- 1
view_images()

view_image_by_name_and_subfolder("20230725_Tree_Shade", 
                                 "Thermal_masked_crop_D42_W_Tree_Shade")


# Create processsing examples and save images
output_dir <- "./graphics/Thermal_Image_Processing_Demo/"  

subfolder <- "20230818_Tree_Sun"
filename_no_ext <- "Thermal_masked_crop_B31_W_Tree_Sun"
filename <- paste0(filename_no_ext, ".tif")

# path to rgb file
out_file <- file.path(output_dir, 
                      paste0("20230818", filename_no_ext, "_plot.png"))

png(out_file, width = 1800, height = 600, res = 150)

# Plot directly
view_image_by_name_and_subfolder(subfolder, filename_no_ext)
dev.off()

