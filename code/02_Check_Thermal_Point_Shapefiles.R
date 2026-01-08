# Check point data from shapefiles created as input for tree segmentation 

library(terra)
library(sf)
library(raster)
library(viridis)

source("./Code/Functions.R")

path = "D:/Bewaesserung_Forstkulturen/Daten/Gewaechshaus/Gewaechshaus_Trockenstress2023/Thermal/GH_Analyse_Thermal/"

# select shapefile-folder and corresponding thermal- and rgb-image-folder
all_dirs <- list.dirs(path, full.names = TRUE, recursive = FALSE)

shp_folderlist <- all_dirs[grep("*_Shapefiles", all_dirs)]
print(shp_folderlist)
shp_folder <- shp_folderlist[[1]]

image_folderlist <- all_dirs[grep("*_Thermal_sorted", all_dirs)]
print(image_folderlist)
image_folder <- image_folderlist[1]
imagesub_folder <- paste0(image_folder, "/Processed_new/")
imagesub_folder_files <- list.files(imagesub_folder, full.names=TRUE)

# select shapefile and corresponding thermal- and rgb- image 

#list all shapefiles in the folder
shps <- list.files(shp_folder, pattern ="*.shp", full.names=TRUE)

# Create a vector of colors for each category of the points 
# (in shp_shift$Position)
category_colors <- c("Tree_SHADE" = "darkred", "Tree_Shade" = "darkred",
                     "Black_SHADE" = "black", "Black_Shade" = "black",
                     "Silver_SHADE" = "darkblue", "Silver_Shade" = "darkblue",
                     "Green_Dry_SHADE" ="orange","Green_Dry_Shade" ="orange",
                     "Green_Wet_SHADE" = "darkgreen",  
                     "Green_Wet_Shade" = "darkgreen",
                     "Tree_SUN" = "pink",  "Tree_Sun" = "pink",
                     "Black_SUN" = "grey", "Black_Sun" = "grey",
                     "Silver_SUN" = "royalblue", "Silver_Sun" = "royalblue",
                     "Green_Dry_SUN" = "yellow","Green_Dry_Sun" = "yellow",
                     "Green_Wet_SUN" = "lightgreen", 
                     "Green_Wet_Sun" = "lightgreen")


# check points in shapefiles
for (i in 1: length(shps)){
  point_check(i)
  print(i)
}


