# Function for splitting dataframe into several dataframes, so one dataframe
# per tree species and treatment is created
tst_df_creation <- function(df_to_split, df_name, treespecies, treatment) {
  selected_columns <- grep(paste0("^", treespecies, ".*", treatment, 
                                  "(_original)?$"), colnames(df_to_split), 
                           value = TRUE)
  dataframe_name <- df_to_split[, c("TIME", selected_columns)]
}

# Function for creating a dataframe per tree species and treatment with 
# dendrometers that have missing data
select_dendro_with_NAs <- function(x) {
  misdat_names_with_suffix <- paste0(dendro_misdat_names, "_original")
  combined_misdat_names <- c(dendro_misdat_names, misdat_names_with_suffix)
  
  if (any(combined_misdat_names %in% colnames(as.data.frame(x)))) {
    x_cleaned <- subset(as.data.frame(x), 
                        select = colnames(as.data.frame(x)) %in% combined_misdat_names)
    return(x_cleaned)
  }
}

# Function for creating a dataframe per tree species and treatment without 
# dendrometers that have missing data
remove_dendro_with_NAs <- function(x, mis_dendro_names) {
  mis_dendro_names_with_suffix <- paste0(mis_dendro_names, "_original")
  combined_mis_dendro_names <- c(mis_dendro_names, mis_dendro_names_with_suffix)
  
  if (any(combined_mis_dendro_names %in% colnames(as.data.frame(x)))) {
    x_cleaned <- subset(as.data.frame(x), 
                        select = !colnames(as.data.frame(x)) %in% combined_mis_dendro_names)
    return(x_cleaned)
  } else {
    return(x)
  }
}
# Function for interpolating missing data with network.interpolation and 
# correcting jumps with dendRoAnalyst (Aryal et al. 2020), returning a dataframe
# with the interpolated and jump-corrected values
process_dendrometer <- function(df1, df2) {
  dendro_network <- network.interpolation(df1, df2, niMethod = "proportional")
  dendro_network_j <- jump.locator(dendro_network, v=30)
  return(cbind(df1, dendro_network, dendro_network_j))
}

# Function for combining the dendRoAnalyst-processed dendrometerdata into one 
# dataframe
create_process_dendro_df <- function(df_column, df){
  dendro_processed <- process_dendrometer(df_column, df)
  col_amount <- ncol(as.data.frame(df_column))
  # remove 3times occurring ts -column
  combined_dendro <- dendro_processed[,c(1:col_amount,
                                         (col_amount+2):(2*col_amount),
                                         (2*col_amount+2):(2*col_amount+2+
                                                             col_amount-2))]
  colnames(combined_dendro) <- c(colnames(combined_dendro[,1:col_amount]),
                                 paste0(colnames(combined_dendro
                                                 [,2:col_amount]), "_int"),
                                 paste0(colnames(combined_dendro
                                                 [,2:col_amount]), "_jcor"))
  return(combined_dendro)
  
}

# Function for correcting one dendrometer timeseries with missing data  
process_single_dendro <- function(df_NA_dendro, corrected_dendro) {
  if (!is.null(df_NA_dendro)) {
    
    # select corrected dendrometer for network interpolation
    df_corrected <- corrected_dendro
    colnames(df_corrected)[1] <- "ts"
    df_na_dendro <- cbind(df_corrected$ts, df_NA_dendro)
    colnames(df_na_dendro)[1] <- "ts"
    df_na_dendro_corrected <- create_process_dendro_df(df_na_dendro, df_corrected)
    
    # change column names
    if (ncol(df_na_dendro_corrected) == 4) {
      colnames(df_na_dendro_corrected) <- c("TIME",
                                            colnames(df_na_dendro_corrected)[2],
                                            paste0(colnames(
                                              df_na_dendro_corrected)[2], 
                                              colnames(df_na_dendro_corrected)[3]),
                                            paste0(colnames(
                                              df_na_dendro_corrected)[2], 
                                              colnames(df_na_dendro_corrected)[4]))
    } else {
      colnames(df_na_dendro_corrected)[1] <- "TIME"
    }
    
    corrected_dendro <- left_join(corrected_dendro, 
                                  df_na_dendro_corrected, by = "TIME")
    
  }
  # Return dataframe with corrected dendrometer
  return(corrected_dendro)
}

# Function for plotting treenetproc corrected, dendRoAnalyst-interpolated and 
# -jump-corrected dendrometer data in comparison per tree.id
plot_dendrometer_comparison <- function(data, tree_col, png_path) {
  # check if a column with dendRoAnalyst-jump-corrected data exists
  jcor_column_exists <- paste0(tree_col, "_jcor") %in% colnames(data)
  
  # create base plot with treenetproc corrected data
  p <- ggplot(data) +
    geom_line(aes(x = Datetime, y = .data[[tree_col]], color = "orig"), 
              linewidth = 0.5) +
    labs(
      title = paste0(tree_col,
                     " Comparison of original
       dendRoAnalyst (network interpolated & jump 30 corrected) processed TSD"),
      x = "Date",
      y = "Tree Stem Diameter [µm]",
      color = ""
    ) +
    scale_color_manual(
      values = c("orig" = "black", 
                 "DA" = "blue"),
      labels = c("orig" = "original", 
                 "DA" = "dendRoAnalyst corrected")
    ) +
    theme_bw() +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),  
          axis.title = element_text(size = 10),  
          axis.text = element_text(size = 7),
          legend.text = element_text(size = 10))
  
  # add dendRoAnalyst-interpolated and -jump-corrected data if it exists
  if (jcor_column_exists) {
    p <- p +
      geom_line(aes(x = Datetime, y = .data[[paste0(tree_col, "_jcor")]], 
                    color = "DA"), 
                linewidth = 0.5) +
      geom_line(aes(x = Datetime, y = .data[[tree_col]], color = "orig"), 
                linewidth = 0.5)
    
  }
  ggsave(plot = p, path = png_path, 
         filename = paste0(tree_col, "_TSD_orig_DA_comparison.png"),
         width = 2560, height = 1489, units = "px")
}

# function for calculating measures like TWD and GRO from dendRoAnalyst and
# create plot with TWD and GRO per Tree.ID
calc_plot_DA <- function (i, df, year, plotpath){
  Tree.ID <- colnames(df)[i]
  Treatment <- gsub(".*([DW]).*", "\\1", Tree.ID)
  Tree.Species <- substr(Tree.ID, 1, 1)
  TWD <- phase.zg(df, TreeNum = i-1)
  ZG_phase <- TWD[[2]]
  ZG_phase$Tree.ID <- Tree.ID
  ZG_phase$Treatment <- Treatment
  ZG_phase$Tree.Species <- Tree.Species
  ZG_phase_DF <<- rbind(ZG_phase_DF, ZG_phase)
  png(paste0(plotpath, Tree.ID, "_TWD.png"), width = 2000, height = 1000)
  p <- plot_ZG_output(ZG_output=TWD,DOY=c(first_day_of_year, last_day_of_year),
                      Year=year)
  print(p)
  dev.off()
}


# Function for cacting the max. log(TWD) per tree species and creating a 
# dataframe with the results where also the Tree.ID and date is listed on which
# it occured
get_max_TWD_log_info <- function(species) {
  max_row <- ZG_phase_DF %>%
    filter(Tree.Species == species) %>%
    filter(TWD_log == max(TWD_log, na.rm = TRUE)) %>%
    slice(1)  
  species2 <- ifelse(max_row$Tree.Species == "Fagus sylvatica", "B",
                     ifelse(max_row$Tree.Species == "Quercus robur", "E", "D"))
  data.frame(
    variable = paste0("max_TWD_", species2, "_log"),
    max_TWD_value = max_row$TWD_log,
    max_TWD_date = max_row$TIME,
    max_TWD_Tree.ID = max_row$Tree.ID
  )
}

# Function for normalizing log(TWD) to max log(TWD) per tree species
normalize_TWD <- function(species) {
  species2 <- ifelse(species == "Fagus sylvatica", "B",
                            ifelse(species == "Quercus robur", "E", "D"))
  
  # extract max log(TWD) according to tree species
  max_value_spec <- max_data_log %>%
    filter(str_detect(variable, paste0("_", species2, "_"))) %>%
    pull(max_TWD_value)
  
  # Normalize log(TWD) to max. log(TWD) per tree species
  ZG_phase_DF %>%
    filter(Tree.Species == species) %>%
    mutate(rel_TWD_log = TWD_log / max_value_spec)
}

# function for checking the manually created points in the shapefiles used 
# for the tree segmentation and the corresponding thermal and rgb images
point_check <- function(x){
  #extract Tree-ID +Treatment from filename
  shape <- shps[[x]]
  Tree_ID_Treat <- gsub(".*0200_(.*)_Points.shp", "\\1", shape)
  
  #read in shapefile with points
  shp <- vect(shape)
  print(shp)
  print(shp$id)
  
  # read in thermal und rgb image
  therm_img <- rast(grep(paste0(Tree_ID_Treat, "_thermal.jpg"), 
                         imagesub_folder_files, value = TRUE))
  
  # flip image, otherwise it will be displayed the wrong way round
  therm_img <- flip(therm_img, direction = "vertical")
  
  print(therm_img)
  
  rgb_img <- rast(grep(paste0(Tree_ID_Treat, "_thermal_rgb_transformed"), 
                       imagesub_folder_files, value = TRUE))
  rgb_img <- flip(rgb_img, direction = "vertical")
  print (rgb_img)
  
  # extract image path
  therm_path <- grep(paste0(Tree_ID_Treat, "_thermal.jpg"), 
                     imagesub_folder_files, value = TRUE)
  rgb_path <- grep(paste0(Tree_ID_Treat, "_thermal_rgb_transformed"), 
                   imagesub_folder_files, value = TRUE)
  
  
  # set coordinate system of thermal and rgb image to the one of the shapefile
  crs(therm_img) <- crs(shp) 
  crs(rgb_img) <-  crs(shp)
  
  # shift point shapefile so it overlays with the thermal- and rgb-image
  shp_shift <- shift(shp, dx= 0, dy= 640)

   # Create a new plot with two side-by-side panels
  windows()
  par(mfrow = c(1, 2))
  
  # Plot the thermal image
  plot(therm_img, main = paste0(Tree_ID_Treat, " Thermal Image"), 
       legend = TRUE, col = inferno(256))

  # Plot the shapefile points with colors based on categories overlaying the 
  # thermal image
  points(shp_shift, pch = 16, col = category_colors[shp_shift$Position])
  
  # Plot the rgb image
  plot(rgb_img, main = paste0(Tree_ID_Treat, " RGB Image"), legend = FALSE)
  
  # Plot the shapefile points with colors based on categories overlaying the 
  # RGB-image
  points(shp_shift, pch = 16, col = category_colors[shp_shift$Position])
  
  # Wait for mouseclick
  click <- locator(1)
  
  # check, if there has been a mouse click
  if (!is.null(click)) {
    # Close Plot window
    dev.off()
    x <- x+1
  }
}

# function for creating and saving rasters (radiance und temperatur)
thermal_raw_rasters <- function(image_folder, crs_template){
  image_files <- list.files(image_folder, pattern = "*.jpg", full.names = TRUE)
  
  for (x in image_files) {
    # Extract metainfo
    imgname <- gsub(".*/(ThermalCamera\\d+-\\d+-\\d+_\\d+-\\d+-\\d+\\+\\d+_\\w+).*", "\\1", x)
    date_val <- gsub(".*ThermalCamera(\\d{4}-\\d{2}-\\d{2}).*", "\\1", imgname)
    date_code <- gsub("-", "", date_val)
    tree_id <- gsub(".*ThermalCamera\\d{4}-\\d{2}-\\d{2}_\\d{2}-\\d{2}-\\d{2}\\+\\d{4}_(\\w+).*", "\\1", imgname)
    
    # check entry in empty_files
    empty_check <- empty_files %>%
      mutate(Date_fmt = format(as.Date(Date, format = "%d.%m.%Y"), "%Y%m%d")) %>%
      filter(Tree.ID == tree_id & Date_fmt == date_code)
    
    if (nrow(empty_check) > 0) {
      message(paste("Skip image for empty entry:", tree_id, date_val))
      next
    }
    rawimg <- readflirJPG(x, exiftoolpath = "installed")
    cams <- flirsettings(x, exiftoolpath = "installed", camvals = "")
    
    # Camera parameter
    ObjectEmissivity <- 0.98
    PlanckR1 <- cams$Info$PlanckR1
    PlanckB <- cams$Info$PlanckB
    PlanckF <- cams$Info$PlanckF
    PlanckO <- cams$Info$PlanckO
    PlanckR2 <- cams$Info$PlanckR2
    ATA1 <- cams$Info$AtmosphericTransAlpha1
    ATA2 <- cams$Info$AtmosphericTransAlpha2
    ATB1 <- cams$Info$AtmosphericTransBeta1
    ATB2 <- cams$Info$AtmosphericTransBeta2
    ATX <- cams$Info$AtmosphericTransX
    OD <- cams$Info$ObjectDistance
    ReflT <- cams$Info$ReflectedApparentTemperature
    AtmosT <- cams$Info$AtmosphericTemperature
    IRWinT <- cams$Info$IRWindowTemperature
    IRWinTran <- cams$Info$IRWindowTransmission
    RH <- cams$Info$RelativeHumidity
    
    h <- cams$Info$RawThermalImageHeight
    w <- cams$Info$RawThermalImageWidth
    
    temperature <- raw2temp(rawimg, ObjectEmissivity, OD, ReflT, AtmosT, IRWinT,
                            IRWinTran, RH, PlanckR1, PlanckB, PlanckF, PlanckO,
                            PlanckR2, ATA1, ATA2, ATB1, ATB2, ATX)
    
    raw_raster <- rast(rawimg, crs = crs_template, extent = c(0, 480, -640, 0))
    temp_raster <- rast(temperature, crs = crs_template, extent = c(0, 480, -640, 0))
    
    date_code <- gsub(".*(\\d{8}).*", "\\1", x)
    raw_folder <- file.path(path, paste0(date_code, "_RawRaster"))
    temp_folder <- file.path(path, paste0(date_code, "_TemperaturRaster"))
    if (!dir.exists(raw_folder)) dir.create(raw_folder)
    if (!dir.exists(temp_folder)) dir.create(temp_folder)
    
    imgname <- gsub(".*/(ThermalCamera\\d+-\\d+-\\d+_\\d+-\\d+-\\d+\\+\\d+_\\w+).*", "\\1", x)
    writeRaster(raw_raster, file.path(raw_folder, paste0(imgname, "_raw.tif")), overwrite = TRUE)
    writeRaster(temp_raster, file.path(temp_folder, paste0(imgname, "_temperature.tif")), overwrite = TRUE)
  }}

# function for extracting temperature and radiance values from corresponding 
# raster files for each point in a shapefile and returns a combined data frame 
# with these values and related metadata (e.g., image name, date, position)
extract_temperatures <- function(shape, raws, temps) {
  # Check if data exists
  if (!file.exists(shape)) {
    message(paste("Shapefile does not exist:", shape))
    return(NULL)
  }
  
  # load shapefile 
  shp <- tryCatch({
    vect(shape)
  }, error = function(e) {
    message(paste("Error reading shapefile:", shape))
    return(NULL)
  })
  
  # Check for empty shapefiles 
  if (nrow(shp) == 0) {
    message(paste("Shapefile does not contain points:", shape))
    
    # Extract meta-Info from filename
    imgname <- gsub(".*(ThermalCamera\\d{4}-\\d{2}-\\d{2}_\\d{2}-\\d{2}-\\d{2}\\+\\d{4}_\\w+).*", "\\1", shape)
    date_val <- gsub(".*ThermalCamera(\\d{4}-\\d{2}-\\d{2}).*", "\\1", imgname)
    tree_id <- gsub(".*ThermalCamera\\d{4}-\\d{2}-\\d{2}_\\d{2}-\\d{2}-\\d{2}\\+\\d{4}_(\\w+)_([WD])_Points.*", "\\1", shape)
    
    # Empty dataframe with NA for temperatur and radiance
    empty_df <- data.frame(
      ID = NA,
      Temperature = NA,
      Radiances_raw = NA,
      Position = NA,
      imgname = imgname,
      Date = date_val,
      Tree.ID = tree_id
    )
    return(empty_df)
  }
  
  Tree_ID_Treat <- gsub(".*0200_(.*)_Points.shp", "\\1", shape, ignore.case = TRUE)
  
  # Load matching raster data (raw values/radiances + temperatures)
  raw_file <- grep(Tree_ID_Treat, raws, value = TRUE)
  temp_file <- grep(Tree_ID_Treat, temps, value = TRUE)
  
  if (length(raw_file) == 0 | length(temp_file) == 0) {
    message(paste("No matching raster for:", Tree_ID_Treat))
    return(NULL)
  }
  
  raw_ras <- tryCatch({
    rast(raw_file)
  }, error = function(e) {
    message(paste("Error reading raw raster:", raw_file))
    return(NULL)
  })
  
  temp_ras <- tryCatch({
    rast(temp_file)
  }, error = function(e) {
    message(paste("Error reading temperature raster:", temp_file))
    return(NULL)
  })
  
  # Extraction: temperatur
  extracted_data_temp <- tryCatch({
    extract(temp_ras, shp)
  }, error = function(e) {
    message("Error extracting temperature")
    return(NULL)
  })
  
  if (is.null(extracted_data_temp) || nrow(extracted_data_temp) == 0) {
    message("Extraction from temperatur raster does not result in data")
    return(NULL)
  }
  
  colnames(extracted_data_temp)[2] <- "Temperature"
  extracted_data_temp$Position <- shp$Position
  
  # Extraction: Rawdata
  extracted_data_rawimg <- tryCatch({
    extract(raw_ras, shp)
  }, error = function(e) {
    message("Error during radiance extraction")
    return(NULL)
  })
  
  if (is.null(extracted_data_rawimg) || nrow(extracted_data_rawimg) == 0) {
    message("Extraction from raw raster does not result in data")
    return(NULL)
  }
  
  colnames(extracted_data_rawimg)[2] <- "Radiances_raw"
  
  # Join and add additional info
  extracted_all <- left_join(extracted_data_temp, extracted_data_rawimg, by ="ID")
  imgname <- gsub(
    ".*(ThermalCamera\\d{4}-\\d{2}-\\d{2}_\\d{2}-\\d{2}-\\d{2}\\+\\d{4}_\\w+).*",
    "\\1", sources(temp_ras))
  extracted_all$imgname <- imgname
  extracted_all$Date <- gsub(
    ".*ThermalCamera(\\d{4}-\\d{2}-\\d{2}).*", "\\1", imgname)
  extracted_all$Tree.ID <- gsub(
    ".*ThermalCamera\\d{4}-\\d{2}-\\d{2}_\\d{2}-\\d{2}-\\d{2}\\+\\d{4}_(\\w+)_([WD])_temperature.*",
    "\\1", imgname)
  return(extracted_all)
}

# function for date extraction
extract_date <- function(folder_name) {
  m <- regmatches(folder_name, regexpr("\\d{8}", folder_name))
  if (length(m) == 1) return(m) else return(NA)
}

# function to extract Tree.ID from filename 
# (between Timestamp and "_thermal_rgb_transformed")
extract_tree_id <- function(filename) {
  pattern <- "_(\\w+_\\w+)_thermal_rgb_transformed"
  m <- regmatches(filename, regexec(pattern, filename))
  if (length(m[[1]]) >= 2) {
    return(m[[1]][2])
  } else {
    return(NA)
  }
}

# function for extracting the Tree.ID from shapefilename
extract_treatment_id <- function(tree_id_full) {
  str_match(tree_id_full, "_([A-Z0-9]+_[A-Z])_Points\\.shp$")[,2]
}

# function for eroding one pixel at the object edge
erode_one_pixel <- function(rast_mask) {
  # Mask for NA-cells
  na_mask <- is.na(rast_mask)
  
  # Neighbouring cells of the NA cells (focal with 3x3 window, max)
  na_expanded <- focal(na_mask, w=matrix(1,3,3), fun=max, na.policy="omit",
                       fill=FALSE)
  
  # Set all cells neihgbouring an NA cell to NA
  rast_mask[na_expanded == 1] <- NA
  return(rast_mask)
}

# function for converting numeric values into RGB color values based on a 
# specified color palette and value range, mapping each value to its 
# corresponding color and returning the result as an RGB matrix
vals_to_rgb <- function(vals, zlim, palette) {
  n_colors <- length(palette)
  idx <- round((vals - zlim[1]) / diff(zlim) * (n_colors - 1)) + 1
  idx[idx < 1] <- 1
  idx[idx > n_colors] <- n_colors
  cols <- rep(NA_character_, length(vals))
  valid <- !is.na(vals)
  cols[valid] <- palette[idx[valid]]
  cols[is.na(cols)] <- "#FFFFFF"  # NA is white
  rgb_mat <- t(col2rgb(cols))
  return(rgb_mat)
}

# function for visually comparing the k-means clustered thermal images to their
# original ones
view_images <- function() {
  col_map <- magma(100)

  for(i in seq(from = start_index, to = length(input_files))) {
    orig_path <- input_files[i]
    print(i)
    rel_path <- gsub(paste0("^", input_root, "/?"), "", orig_path)
    kmeans_path <- file.path(target_root, rel_path)

    if (!file.exists(kmeans_path)) {
      message("Missing KMeans-image: ", rel_path)
      next
    }

    img_orig <- rast(orig_path)
    img_kmeans <- rast(kmeans_path)

    zlim <- range(c(values(img_orig), values(img_kmeans)), na.rm=TRUE)

    vals_orig <- values(img_orig)
    vals_kmeans <- values(img_kmeans)

    rgb_orig <- vals_to_rgb(vals_orig, zlim, col_map)
    rgb_kmeans <- vals_to_rgb(vals_kmeans, zlim, col_map)

    r_orig_rgb <- rast(nlyrs=3, nrow=nrow(img_orig), ncol=ncol(img_orig),
                       extent=ext(img_orig), crs=crs(img_orig))
    r_kmeans_rgb <- rast(nlyrs=3, nrow=nrow(img_kmeans), ncol=ncol(img_kmeans),
                         extent=ext(img_kmeans), crs=crs(img_kmeans))

    values(r_orig_rgb[[1]]) <- rgb_orig[,1]
    values(r_orig_rgb[[2]]) <- rgb_orig[,2]
    values(r_orig_rgb[[3]]) <- rgb_orig[,3]

    values(r_kmeans_rgb[[1]]) <- rgb_kmeans[,1]
    values(r_kmeans_rgb[[2]]) <- rgb_kmeans[,2]
    values(r_kmeans_rgb[[3]]) <- rgb_kmeans[,3]

    layout(matrix(c(1,2,3), 1, 3), widths=c(1,1,0.3))
    par(mar=c(3,3,4,2))

    plotRGB(r_orig_rgb, r=1, g=2, b=3, axes=FALSE)
    title(main=paste0("Original Thermal Image\n", basename(orig_path)), line=2)

    plotRGB(r_kmeans_rgb, r=1, g=2, b=3, axes=FALSE)
    title(main=paste0("Processed Thermal Image\n", basename(kmeans_path)), line=2)

    par(mar=c(3,1,4,4))
    n_col <- length(col_map)
    y_vals <- seq(zlim[1], zlim[2], length.out = n_col + 1)
    z_mat <- matrix(seq(zlim[1], zlim[2], length.out = n_col), nrow = 1, ncol = n_col)

    image(x = c(0,1), y = y_vals, z = z_mat,
          col = col_map, axes = FALSE)
    axis(4, at = pretty(y_vals), labels = pretty(y_vals), las = 1, cex.axis = 0.8)
    title("Temperature", line = 2.5)



    readline(prompt="Press [Enter] for next image...")
  }
  message("Completed.")
}


## view selected image
view_image_by_name_and_subfolder <- function(subfolder, 
                                             filename_no_ext, 
                                             output_dir = NULL) {
  col_map <- viridis::magma(100)
  filename <- paste0(filename_no_ext, ".tif")
  
  # --- Original Thermal ---
  match_path <- input_files[
    grepl(subfolder, input_files) & grepl(filename, basename(input_files))
  ]
  if (length(match_path) == 0) stop("Image not found: ", filename, 
                                    " in folder ", subfolder)
  orig_path <- match_path[1]
  
  # --- KMeans Thermal ---
  rel_path <- gsub(paste0("^", input_root, "/?"), "", orig_path)
  kmeans_path <- file.path(input_root_kmeans, rel_path)
  if (!file.exists(kmeans_path)) stop("Missing KMeans-image: ", kmeans_path)
  
  # --- Processed Thermal ---
  proc_path <- file.path(target_root, rel_path)
  if (!file.exists(proc_path)) stop("Missing processed image: ", proc_path)
  
  # load raster
  img_orig   <- rast(orig_path)
  img_kmeans <- rast(kmeans_path)
  img_proc   <- rast(proc_path)
  
  # common color scale
  zlim <- range(c(values(img_orig), values(img_kmeans), 
                  values(img_proc)), na.rm=TRUE)
  
  # change numeric values to rgb values
  rgb_orig   <- vals_to_rgb(values(img_orig),   zlim, col_map)
  rgb_kmeans <- vals_to_rgb(values(img_kmeans), zlim, col_map)
  rgb_proc   <- vals_to_rgb(values(img_proc),   zlim, col_map)
  
  # Create RGB-Raster 
  make_rgb_rast <- function(img, rgb_vals) {
    r <- rast(nlyrs=3, nrow=nrow(img), ncol=ncol(img), extent=ext(img), crs=crs(img))
    values(r[[1]]) <- rgb_vals[,1]
    values(r[[2]]) <- rgb_vals[,2]
    values(r[[3]]) <- rgb_vals[,3]
    return(r)
  }
  
  r_orig_rgb   <- make_rgb_rast(img_orig,   rgb_orig)
  r_kmeans_rgb <- make_rgb_rast(img_kmeans, rgb_kmeans)
  r_proc_rgb   <- make_rgb_rast(img_proc,   rgb_proc)
  
  # Extract Tree.ID, object, date 
  folder_name    <- basename(dirname(orig_path))
  date_str       <- gsub("^(\\d{8}).*$", "\\1", folder_name)
  date_formatted <- format(as.Date(date_str, format="%Y%m%d"), "%d.%m.%Y")
  file_base      <- tools::file_path_sans_ext(basename(orig_path))
  tree_id        <- gsub(".*_((B|D)[0-9]+_[DW]).*", "\\1", file_base)
  object         <- gsub(".*_(Tree_Shade|Tree_Sun|Green_Dry_Shade|Green_Dry_Sun)$", 
                         "\\1", file_base)
  object         <- gsub("_", " ", object)
  title_text     <- paste(tree_id, object, "\n", date_formatted)
  
  # Layout: Original | KMeans | Processed | Farbskala
  layout(matrix(c(1,2,3,4), 1, 4), widths = c(1,1,1,0.2))
  title_cex  <- 1.7
  title_line <- 2.5
  top_margin <- 6
  
  # --- Panel 1: Original ---
  par(mar = c(3, 3, top_margin, 2))
  plotRGB(r_orig_rgb, r=1, g=2, b=3, axes=FALSE)
  title(main = "Thermal Image", line=title_line, cex.main=title_cex, font.main=2)
  title(sub  = title_text, line=0.8, cex.sub=1.3, font.sub=2)
  
  # --- Panel 2: KMeans ---
  par(mar = c(3, 3, top_margin, 2))
  plotRGB(r_kmeans_rgb, r=1, g=2, b=3, axes=FALSE)
  title(main = "KMeans Thermal", line=title_line, cex.main=title_cex, font.main=2)
  title(sub  = title_text, line=0.8, cex.sub=1.3, font.sub=2)
  
  # --- Panel 3: Processed ---
  par(mar = c(3, 3, top_margin, 2))
  plotRGB(r_proc_rgb, r=1, g=2, b=3, axes=FALSE)
  title(main = "Processed Thermal", line=title_line, cex.main=title_cex, font.main=2)
  title(sub  = title_text, line=0.8, cex.sub=1.3, font.sub=2)
  
  # --- Panel 4: Color scale ---
  par(mar = c(3, 1, top_margin, 4))
  n_col  <- length(col_map)
  y_vals <- seq(zlim[1], zlim[2], length.out = n_col + 1)
  z_mat  <- matrix(seq(zlim[1], zlim[2], length.out = n_col), nrow=1, ncol=n_col)
  
  image(x=c(0,1), y=y_vals, z=z_mat, col=col_map, axes=FALSE)
  axis(4, at=pretty(y_vals), labels=pretty(y_vals), las=1, cex.axis=1.3)
  mtext("Temperature [°C]", side=3, line=title_line, cex=title_cex-0.6, font=2)
}

# Helperfunction for performance metrics from crossvalidation
extract_metrics <- function(preds, true) {
  cm <- caret::confusionMatrix(preds, true)
  sens <- cm$byClass["Sensitivity"]   # Recall for positive class
  spec <- cm$byClass["Specificity"]   # Recall for negative class
  bal_acc <- mean(c(sens, spec))      # Balanced Accuracy
  data.frame(
    Accuracy  = cm$overall["Accuracy"],
    Kappa     = cm$overall["Kappa"],
    Precision = cm$byClass["Pos Pred Value"],
    Recall    = cm$byClass["Sensitivity"],
    F1        = cm$byClass["F1"],
    Balanced_Accuracy = bal_acc
  )
}

# Function for formatting p-values
format_pval <- function(p) {
  if (is.na(p)) return(NA)
  if (p < 0.001) return("< 0.001")
  if (p < 0.01) return("< 0.01")
  if (p < 0.05) return("< 0.05")
  return(round(p, 3))
}

# helper function for boxplots representing differences in CWSI, SWC and VPD
# according to TWDmin classes
make_boxplot <- function(df, yvar, ylab, show_x = FALSE) {
  ggplot(df, aes(x = Tree.Species, y = .data[[yvar]], fill = class)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.size = 2,
                 color = "black", width = 0.5,
                 show.legend = show_x) +   # only lowest plot has legend
    
    stat_compare_means(
      aes(group = class),
      method = "wilcox.test",
      label = "p.signif",
      label.y = max(df[[yvar]], na.rm = TRUE) * 0.95,
      hide.ns = TRUE,
      size = 6
    ) +
    
    labs(
      x = if(show_x) "Tree species" else NULL,
      y = ylab,
      fill = if(show_x) expression(TWD[min]) else NULL
    ) +
    
    scale_fill_manual(values = cols) + 
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = if(show_x)
        element_text(angle = 0, hjust = 0.5, size = 13) else element_blank(),
      axis.ticks.x = if(show_x) element_line() else element_blank(),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold"),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      legend.key.height = unit(0.8, "cm"),
      legend.spacing.y = unit(0.3, "cm"),
      plot.margin = margin(t = 15, r = 10, b = 15, l = 10)
    )
}
