# Create temperature rasters from the raw values for every image and
# extract temperature values from the selected points in the shapefiles 


library(terra)
library(viridis)
library(Thermimage)
library(fields)
library(dplyr)
library(readr)
library(stringr)


source("./Code/Functions.R")

path = "D:/Bewaesserung_Forstkulturen/Daten/Gewaechshaus/Gewaechshaus_Trockenstress2023/Thermal/GH_Analyse_Thermal/"

image_folderlist <- list.dirs(path, full.names = TRUE, recursive = FALSE)
image_folderlist <- image_folderlist[grep("Thermal_sorted", image_folderlist)]

# info with dates and tree ids when an empty image shall be created, since there
# was no thermal image taken
empty_files <- read_csv(
  "./data/Thermal/Empty_Images.csv")
empty_files$Date <- as.POSIXct(empty_files$Date, format = "%d.%m.%Y")

# Iterate through all folders
for (image_folder in image_folderlist) {
  print(image_folder)
  date_code <- gsub(".*(\\d{8}).*", "\\1", image_folder)
  print(paste("Proceed:", date_code))
  
  # read shapefile
  shp_folder <- file.path(path, paste0(date_code, "_Shapefiles"))
  print(shp_folder)
  shps <- list.files(shp_folder, pattern = ".*[Pp]oints.shp$", full.names = TRUE)
  
  if (length(shps) == 0) {
    warning(paste("No shapefiles for", date_code, "found – skipped."))
    next
  }
  
  # Load example shapefile to extract cordinate reference system
  shp_example <- vect(shps[1])
  
  # Create raster
  thermal_raw_rasters(image_folder, crs(shp_example))
  
  # Load raster
  raw_folder <- file.path(path, paste0(date_code, "_RawRaster"))
  temp_folder <- file.path(path, paste0(date_code, "_TemperaturRaster"))
  print(raw_folder)
  print(temp_folder)
  
  raws <- list.files(raw_folder, pattern = "*.tif", full.names = TRUE)
  temps <- list.files(temp_folder, pattern = "*.tif", full.names = TRUE)
  
  # Extract values
  result_df <- do.call(bind_rows, 
                       lapply(shps, 
                              function(s) extract_temperatures(s, 
                                                               raws = raws, 
                                                               temps = temps)))
  
  # Correct labeling
  result_df$Position <- sub("SHADE", "Shade", result_df$Position, 
                            ignore.case = TRUE)
  result_df$Position <- sub("SUN", "Sun", result_df$Position, 
                            ignore.case = TRUE)
  
  # add metadata
  metainfos <- read.csv("./data/Tree_Metainfos.csv")[2:7]
  result_df_final <- left_join(result_df, metainfos, by = "Tree.ID")
  
  # Save results as csv
  write.csv(result_df_final, file = paste0("./data/Thermal/tempvalue_tables/",
                                           date_code, 
                                           "_Temperatures+Radiances.csv"))
}

# table entries are empty for the empty shapefiles where do not have data,
# but the temprasters and rawrasters where created and contain some wrong data,
# thus delete them manually


