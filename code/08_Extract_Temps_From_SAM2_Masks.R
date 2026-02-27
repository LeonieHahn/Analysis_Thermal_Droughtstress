# Processes thermal imagery by applying segmentation masks derived from SAM 2
# to isolate tree saplings and reference surfaces in temperature rasters

library(tidyverse)
library(terra)
library(jpeg)
library(stringr)
library(imager)
library(lubridate)
library(doParallel)
library(foreach)
library(doSNOW)

source("./Code/Functions.R")

# calculate in parallel on several cores
num_cores <- 3
cl <- makeCluster(num_cores)
registerDoSNOW(cl)


# Base directories
base_rgb_path <- "./data/images/segmentation_results"
base_thermal_path <- "./data/images/"
output_base <- "./data/images/Cropped_Thermals_1eroded_pixel"

out_dir_allpix <- "./data/images/Cropped_Thermals_all_pixels"
  
# Load table with information which SAM 2 mask should be used
segment_info <- read_csv("./data/Thermal/TreeSpecies_Dates_TreeIDs.csv")
segment_info$Date <- as.Date(as.character(segment_info$Date), format = "%Y%m%d")

# Function for loading the raster with the temperature information, 
# masking the objects (tree sapling & reference surfaces) according to the 
# segmentation masks derived from SAM 2 and savethe masked temperature raster
process_entry <- function(tree_core_id, mask_type, date, channel, 
                          do_allpix = FALSE) {
  if (is.na(channel) || channel == "") return(NULL)
  date_str <- format(date, "%Y%m%d")
  
  # directories
  rgb_dir <- file.path(base_rgb_path, 
                       paste0(date_str, "_Thermal_sorted_", mask_type),
                       channel)
  thermal_dir <- file.path(base_thermal_path, paste0(date_str, "_TemperaturRaster"))
  output_dir <- file.path(output_base, paste0(date_str, "_", mask_type))
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Load RGB and thermal files
  rgb_files <- list.files(rgb_dir, pattern = "\\.jpg$", full.names = TRUE)
  thermal_files <- list.files(thermal_dir, pattern = "\\.tif$", full.names = TRUE)
  
  rgb_file <- rgb_files[str_detect(rgb_files, tree_core_id)]
  thermal_file <- thermal_files[str_detect(thermal_files, tree_core_id)]
  
  if (length(rgb_file) == 0 || length(thermal_file) == 0) {
    cat("File not found for", tree_core_id, mask_type, "\n")
    return(NULL)
  }
  
  cat("Process ", tree_core_id, "–", mask_type, "
      (", date_str, "/", channel, ")\n")
  
  # Create mask
  rgb <- readJPEG(rgb_file[1])
  # remove also darkgrey pixels around objects that are sometimes still included,
  # as the edges are not recognizes too well
  mask_mat <- (rgb[,,1] > 20/255) | (rgb[,,2] > 20/255) | (rgb[,,3] > 20/255)
  nrows <- dim(mask_mat)[1]; ncols <- dim(mask_mat)[2]
  
  rgb_template <- rast(nrows = nrows, ncols = ncols)
  thermal <- rast(thermal_file[1])
  ext(rgb_template) <- ext(thermal)
  crs(rgb_template) <- crs(thermal)
  thermal_resampled <- resample(thermal, rgb_template, method = "bilinear")
  
  mask_rast <- rast(matrix(as.integer(mask_mat), nrow = nrows, ncol = ncols))
  ext(mask_rast) <- ext(thermal_resampled)
  crs(mask_rast) <- crs(thermal_resampled)
  mask_rast[mask_rast < 1] <- NA
  mask_rast[mask_rast >= 1] <- 1
  
  patch_rast <- patches(mask_rast, directions = 8)
  fl <- freq(patch_rast)
  valid_ids <- fl$value[fl$count >= 50]
  
  cleaned_mask <- patch_rast
  cleaned_mask[!patch_rast %in% valid_ids] <- NA
  cleaned_mask[!is.na(cleaned_mask)] <- 1
  
  # crop the thermal image to the SAM2 mask
  
  if (do_allpix) {
    thermal_masked_allpix <- mask(thermal_resampled, cleaned_mask)
    
    out_file_allpix <- file.path(out_dir_allpix,
                                 paste0("Thermal_masked_crop_",
                                        tree_core_id, "_", mask_type, ".tif"))
    writeRaster(thermal_masked_allpix, out_file_allpix, overwrite = TRUE)
  }
  
  
  # Remove pixels at the edge of the mask
  if (mask_type %in% c("Green_Wet_Shade", "Green_Dry_Shade",
                       "Green_Wet_Sun", "Green_Dry_Sun",
                       "Tree_Shade", "Tree_Sun")) {

    # Performe the erosion with 1 pixel (but could also be more pixels),
    # center is kept
    for (i in 1:1) {
      cleaned_mask <- erode_one_pixel(cleaned_mask)
    }
   }
  
  
  thermal_masked <- mask(thermal_resampled, cleaned_mask)
  
  out_file <- file.path(output_dir, paste0("Thermal_masked_crop_", 
                                          tree_core_id, "_", mask_type, ".tif"))
  writeRaster(thermal_masked, out_file, overwrite = TRUE)
}

# Process all 
mask_cols <- c("Green_Dry_Shade", "Green_Dry_Sun", "Green_Wet_Shade", 
               "Green_Wet_Sun", "Tree_Shade", "Tree_Sun")

# prepare all combinations of objects and light conditions, tree.id and date
# for processing
task_list <- list()

for (i in 1:nrow(segment_info)) {
  row <- segment_info[i, ]
  tree_id_full <- row$Tree_ID
  date <- row$Date
  tree_core_id <- extract_treatment_id(tree_id_full)
  
  for (mask_col in mask_cols) {
    channel_raw <- row[[mask_col]]
    if (!is.na(channel_raw) && channel_raw != "") {
      channel <- if (tolower(
        channel_raw) == "combined") "combined" else paste0("channel", 
                                                           channel_raw)
      task_list[[length(task_list) + 1]] <- list(tree_core_id = tree_core_id,
                                                 mask_col = mask_col,
                                                 date = date,
                                                 channel = channel)
    }
  }
}

chunk_size <- 10  # Amount of tasks per chunk
num_chunks <- ceiling(length(task_list) / chunk_size)

# Create process bar (chunk level)
pb_total <- txtProgressBar(min = 0, max = num_chunks, style = 3)

for (chunk in 1:num_chunks) {
  cat("\n Start Chunk", chunk, "from", num_chunks, "\n")
  
  # Select tasks for this junk
  start <- (chunk - 1) * chunk_size + 1
  end <- min(chunk * chunk_size, length(task_list))
  task_subset <- task_list[start:end]
  
  # Processbar for this chunk (task level)
  pb_chunk <- txtProgressBar(min = 0, max = length(task_subset), style = 3)
  progress <- function(n) setTxtProgressBar(pb_chunk, n)
  opts <- list(progress = progress)
  
  foreach(i = seq_along(task_subset),
          .packages = c("terra", "jpeg", "stringr", "lubridate"),
          .options.snow = opts) %dopar% {
            task <- task_subset[[i]]
            tryCatch({
              process_entry(task$tree_core_id, task$mask_col, 
                            task$date, task$channel, do_allpix = FALSE)
            }, error = function(e) {
              cat("Error at:", task$tree_core_id, task$mask_col, "\n")
              NULL
            })
          }
  
  close(pb_chunk)
  setTxtProgressBar(pb_total, chunk)  # Total progress bar
  gc()  # Storage cleaning after every chunk
}

close(pb_total)
stopCluster(cl)
cat("\n Processing done.\n")

for (task in task_list) {
  process_entry(task$tree_core_id,
                task$mask_col,
                task$date,
                task$channel,
                do_allpix = TRUE)
}

cat("\n Processing all pixels done.\n")
