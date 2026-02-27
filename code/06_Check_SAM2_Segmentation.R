# check the segmentation results from SAM 2 for the different objects and 
# light conditions

library(png)
library(grid)
library(gridExtra)

source("./Code/Functions.R")

base_path <- "D:/Bewaesserung_Forstkulturen/Daten/Gewaechshaus/Gewaechshaus_Trockenstress2023/Thermal/GH_Analyse_Thermal/segmentation_results"
all_folders <- list.files(base_path, full.names = TRUE)

all_images <- data.frame(folder = character(),
                         filepath = character(),
                         date = character(),
                         tree_id = character(),
                         stringsAsFactors = FALSE)

# collect metadata from all overview plot files and save in a table
for (folder in all_folders) { 
  files <- list.files(folder, pattern = "\\.jpg_overview.png", full.names = TRUE)
  if (length(files) == 0) next
  
  date <- extract_date(basename(folder))
  
  for (f in files) {
    fname <- basename(f)
    tid <- extract_tree_id(fname)
    if (!is.na(tid)) {
      all_images <- rbind(all_images, data.frame(folder = basename(folder),
                                                 filepath = f,
                                                 date = date,
                                                 tree_id = tid,
                                                 stringsAsFactors = FALSE))
    }}}

dates <- unique(all_images$date)

# plot SAM 2 results of every tree sapling and reference per light condition
# (sunny/shady) for checking which mask shall be used for further analysis
for (current_date in dates) {
  cat("Date:", current_date, "\n")
  subset_date <- subset(all_images, date == current_date)
  tree_ids <- unique(subset_date$tree_id)
  
  for (tid in tree_ids) {
    cat("Tree-ID:", tid, "\n")
    subset_tid <- subset(subset_date, tree_id == tid)
    
    images_list <- list()
    titles <- character()
    
    for (i in seq_len(nrow(subset_tid))) {
      img_file <- subset_tid$filepath[i]
      img <- readPNG(img_file)
      images_list[[length(images_list) + 1]] <- rasterGrob(img, interpolate=TRUE)
      titles <- c(titles, subset_tid$folder[i])
    }
    
    if (length(images_list) > 0) {
      grid.arrange(grobs = images_list, ncol = 2,
                   top = paste0("Date: ", current_date, " | Tree-ID: ", tid),
                   bottom = "Proceed with enter...")
    } else {
      cat("No images for this Tree-ID.\n")
    }
    readline(prompt = "Press enter for next image...")
    grid.newpage()
  }
}

