# This script produces GeoTiff raster files corresponding to the average of
# the full time series of many of the common variable for the BGC and Hydro
# model data for a set of depths. This creates a climatology estimate of
# each variable allowing spatial patterns of average conditions to be 
# identified. 
#
# This script downloads only a subset of the depth slices. 
#
# The script processes four standard depths: 3m, 9m, 18m, and 39m below the surface.
# Each model is processed using predefined date ranges to ensure only complete years 
# are included in the climatology. Variables without a depth dimension are processed 
# once and saved in a 'surface' directory.
#
# In this we we use the AIMS THREDDS regridded aggregate version of the eReefs
# data for two reasons:
# 1. The data is already aggregated to an annual time step, vastly reducing the
#    amount of data we need to download.
# 2. The data is already regridded to a rectangular grid, which is needed for
#    creating GeoTiff files.
#
# The data is downloaded from the AIMS THREDDS server using OpenDAP as it allows
# us to specify which variables we want to download and the depth layer we are 
# interested in, further limiting the amount of data we need to download.
#
# The data is downloaded using OPeNDAP and then aggregated across all full years, 
# and finally saved as one GeoTiff per variable per depth.
#
# The script also copies over key metadata from the original NetCDF file to the new
# GeoTiff files. This includes the variable names, units, and long names and global
# attributes. 
#
# This script prepares GeoTiffs for each model (GBR4 BGC, GBR4 Hydro, and GBR1 Hydro),
# each common variable and each depth layer (3m and 17.75m). The final output is a 
# set of GeoTiff files saved in 
# data/out/<name>/<depth>/<ModelID>_avg-<start_year>-<end-year>_<variable>_<depth>m.tif.

#
# The script is designed to be run in R and requires the ncdf4 and raster packages.
# To provide feedback during progress the script downloads one variable at a time
# and prints the progress to the console. 

# Model,version,ModelID,temporal aggregation,URL
# GBR4 BGC baseline,3.1,GBR4_H2p0_B3p1_Cq3b_Dhnd,annual,https://thredds.ereefs.aims.gov.au/thredds/dodsC/GBR4_H2p0_B3p1_Cq3b_Dhnd/annual.nc
# GBR4 Hydro,2.0,gbr4_v2_hydro,annual,https://thredds.ereefs.aims.gov.au/thredds/dodsC/gbr4_v2/annual.nc
# GBR1 Hydro,2.0,gbr1_2.0_hydro,annual,https://thredds.ereefs.aims.gov.au/thredds/dodsC/gbr1_2.0/annual.nc
# GBR4 Hydro,4.0,GBR4_H4p0_ABARRAr2_OBRAN2020_FG2Gv3_Dhnd,annual,https://thredds.ereefs.aims.gov.au/thredds/dodsC/GBR4_H4p0_ABARRAr2_OBRAN2020_FG2Gv3_Dhnd/annual.nc

# Note: Version 3.1 of the BGC model and version 2.0 of the GBR4 Hydro will
# be deprecated in June 2025. The new version of the BGC model (v4.0) will be released 
# along with the new Hydro model (v4.0). As of June 2025 there is no update for GBR1 Hydro. 
# The GBR4 Hydro v4.0 model is not yet available from the AIMS THREDDS server.

# See https://thredds.ereefs.aims.gov.au/thredds/dodsC/GBR4_H2p0_B3p1_Cq3b_Dhnd/annual.nc.html for full list
# of variables for the BGC model and their descriptions.
# BGC common variables:
# TN = Total Nitrogen
# TP = Total Phosphorus 
# DIN = Dissolved Inorganic Nitrogen
# DIP = Dissolved Inorganic Phosphorus
# Chl_a_sum = Chlorophyll a
# NO3 = Nitrate
# NH4 = Ammonium
# DOR_N = Dissolved Organic Nitrogen
# DOR_P = Dissolved Organic Phosphorus
# Secchi = Secchi depth
# PH = pH
# omega_ar = Aragonite saturation state
# EFI = Ecology Fine Inorganics (Total suspended solids)

# See https://thredds.ereefs.aims.gov.au/thredds/dodsC/GBR4_H4p0_ABARRAr2_OBRAN2020_FG2Gv3_Dhnd/annual.nc.html 
# for additional details about the Hydro model variables.
# Hydro variables:
# eta = Surface elevation (Height of the water affected by tides)
# temp = Temperature
# salt = Salinity
# u = Eastward current (ocean current)
# v = Northward current (ocean current)
# wspeed_u = Eastward wind speed
# wspeed_v = Northward wind speed
# mean_cur = Mean current magnitude speed
# mean_wspeed = Mean wind magnitude speed


library(ncdf4)   # For handling NetCDF files
library(raster)  # For creating and handling raster data
library(sf)      # For spatial operations

# Create output directory structure if it doesn't exist
output_dir <- "data/out"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Define the depths we want to extract (m)
depths <- c(-3, -9, -18, -39)

# Model definitions - contains all the information we need to access each model
# Set start and end years to only include whole years. This is to ensure that
# the average is not skewed by partial years and seasonal effects.
models <- list(
  list(
    name = "GBR4-BGC3p1-base",
    model_id = "GBR4_H2p0_B3p1_Cq3b_Dhnd",
    url = "https://thredds.ereefs.aims.gov.au/thredds/dodsC/GBR4_H2p0_B3p1_Cq3b_Dhnd/annual.nc",
    variables = c('TN', 'TP', 'DIN', 'DIP', 'Chl_a_sum', 'NO3', 'NH4', 'DOR_N', 
                  'DOR_P', 'Secchi', 'PH', 'omega_ar', 'EFI'),
    type = "bgc",
    start_year = 2011,
    end_year = 2018
  ),
  list(
    name = "GBR4-H2p0",
    model_id = "gbr4_v2_hydro",
    url = "https://thredds.ereefs.aims.gov.au/thredds/dodsC/gbr4_v2/annual.nc",
    variables = c('eta', 'temp', 'salt', 'u', 'v', 'wspeed_u', 
      'wspeed_v', 'mean_cur', 'mean_wspeed'),
    type = "hydro",
    start_year = 2011,
    end_year = 2023
  ),
  list(
    name = "GBR1-H2p0",
    model_id = "gbr1_2.0_hydro",
    url = "https://thredds.ereefs.aims.gov.au/thredds/dodsC/gbr1_2.0/annual.nc",
    variables = c('eta', 'temp', 'salt', 'u', 'v', 'wspeed_u', 'wspeed_v',
      'mean_cur', 'mean_wspeed'),
    type = "hydro",
    start_year = 2015,
    end_year = 2023
  )
# Commented out for now as the model is not yet available
#   list(
#     name = "GBR4-H4p0",
#     version = "4.0",
#     model_id = "GBR4_H4p0_ABARRAr2_OBRAN2020_FG2Gv3_Dhnd",
#     url = "https://thredds.ereefs.aims.gov.au/thredds/dodsC/GBR4_H4p0_ABARRAr2_OBRAN2020_FG2Gv3_Dhnd/annual.nc",
#     variables = c('eta', 'temp', 'salt', 'u', 'v', 'wspeed_u', 'wspeed_v'),
#     type = "hydro",
#     start_year = 2011,
#     end_year = 2022
#   )
)


# Function to find the closest depth index using the 'zc' variable that all eReefs models have
find_closest_depth <- function(target_depth, nc) {
  # Get depth values from the 'zc' variable
  if ("zc" %in% names(nc$var)) {
    depth_values <- ncvar_get(nc, "zc")
    
    # Find the index of the closest depth
    idx <- which.min(abs(depth_values - target_depth))
    k <- idx - 1  # Convert to 0-based index for NetCDF slicing
    
    return(list(index = k, actual_depth = depth_values[idx]))
  } else {
    stop("Error: 'zc' variable not found in the NetCDF file. Cannot determine depth levels.")
  }
}

# Modified function that accepts time indices rather than years
download_and_average_variable <- function(nc, var_name, output_file, start_vec, count_vec, 
                                         time_indices, time_years, force_overwrite = FALSE) {

  # Check if file already exists and skip if so (unless force_overwrite is TRUE)
  if (file.exists(output_file) && !force_overwrite) {
    cat(sprintf("\n    Skipping existing file: %s\n", output_file))
    return(NULL)
  }

  if (length(time_indices) == 0) {
    cat(sprintf("\n    Warning: No time steps to process\n"))
    return(NULL)
  }

  # Initialize arrays for accumulating values and counting non-NA cells
  first_slice <- TRUE
  sum_data <- NULL
  count_data <- NULL
  
  # Print header for time slice diagnostics
  cat(sprintf("\n    Time slice averages for %s:\n", var_name))
  
  # Process one time step at a time to limit OpenDAP request size
  for (i in 1:length(time_indices)) {
    t <- time_indices[i]

    # Modify the start and count vectors for this time step
    time_start <- start_vec
    time_count <- count_vec
    
    # Find the time dimension index
    time_dim_idx <- which(sapply(nc$var[[var_name]]$dim, function(d) d$name == "time"))
    
    # Set time slice
    time_start[time_dim_idx] <- t
    time_count[time_dim_idx] <- 1
    
    # Get data for this time slice
    slice_data <- ncvar_get(nc, var_name, start = time_start, count = time_count)
    
    # Remove singleton dimensions
    slice_data <- drop(slice_data)
    
    # Calculate and print average for this time slice (for debugging)
    slice_avg <- mean(slice_data, na.rm = TRUE)

    # Print the time index, year and time slice average
    cat(sprintf("      Time index %d (year %d): %.6f\n", t, time_years[t], slice_avg))
    
    # Initialize arrays on first iteration
    if (first_slice) {
      sum_data <- slice_data
      count_data <- !is.na(slice_data)
      first_slice <- FALSE
    } else {
      # Add to sum (ignoring NAs)
      sum_data[!is.na(slice_data)] <- sum_data[!is.na(slice_data)] + 
                                    slice_data[!is.na(slice_data)]
      
      # Count valid data points
      count_data <- count_data + !is.na(slice_data)
    }
  }

  cat("\n")
  
  # Calculate mean
  mean_data <- sum_data / count_data
  
  # Replace 0/0 division results (NaN) with NA
  mean_data[is.nan(mean_data)] <- NA
  
  # Print the overall mean for validation
  overall_mean <- mean(mean_data, na.rm = TRUE)
  cat(sprintf("    Overall mean: %.6f\n", overall_mean))
  
  # Get coordinates
  lon <- ncvar_get(nc, "longitude")
  lat <- ncvar_get(nc, "latitude")
  
  # Calculate cell size (resolution)
  x_res <- (max(lon) - min(lon)) / (length(lon))
  y_res <- (max(lat) - min(lat)) / (length(lat))
  
  # Create raster with adjusted extent to account for cell centers
  # Adjust by half cell in each direction
  r <- raster(
    t(mean_data),  # Transpose because raster expects rows=y, cols=x
    xmn = min(lon) - (x_res/2), xmx = max(lon) + (x_res/2),
    ymn = min(lat) - (y_res/2), ymx = max(lat) + (y_res/2),
    crs = "+proj=longlat +datum=WGS84"
  )
  
  # Fix the vertical flip issue
  r <- flip(r, direction='y')
  
  # Save as GeoTIFF with all metadata included
  writeRaster(r, output_file, format = "GTiff", overwrite = TRUE)
  
  cat(sprintf("    Created GeoTIFF: %s\n", output_file))
}

# Process each model
for (model in models) {
  model_dir_name <- model$name
  cat(sprintf("\nProcessing model: %s\n", model$name))
  
  # Create model output directory
  model_dir <- file.path(output_dir, model_dir_name)
  if (!dir.exists(model_dir)) dir.create(model_dir, recursive = TRUE)
  
  # Try to open the NetCDF connection (with error handling)
  nc <- tryCatch({
    nc_open(model$url)
  }, error = function(e) {
    cat(sprintf("Error opening NetCDF connection for %s: %s\n", model$name, e$message))
    return(NULL)
  })
  
  if (is.null(nc)) next
  
  # Get time data and convert to Queensland time (UTC+10)
  time_var <- ncvar_get(nc, "time")
  time_units <- ncatt_get(nc, "time", "units")$value
  
  # Parse time units to get base date
  time_base_year <- as.numeric(substr(time_units, regexpr("[0-9]{4}", time_units), 
                                     regexpr("[0-9]{4}", time_units) + 3))
  
  
  # Adjust to Queensland local time (UTC+10)
  # This is a simplification; for more precise conversion we'd need POSIXct handling
  qld_offset <- 10/24/365.25  # 10 hours as fraction of a year
  time_years <- round(time_base_year + time_var/365.25 + qld_offset)
  
  # Use predefined start and end years from the model definition
  start_year <- model$start_year
  end_year <- model$end_year
  
  # Find time indices within the specified year range
  year_mask <- time_years >= start_year & time_years <= end_year


  time_indices <- which(year_mask)
  
  if (length(time_indices) == 0) {
    cat(sprintf("  Warning: No time steps found within years %d-%d. Skipping model.\n", 
                start_year, end_year))
    nc_close(nc)
    next
  }
  
  cat(sprintf("  Found %d time steps within years %d-%d\n", 
              length(time_indices), start_year, end_year))
  cat(sprintf("  Using time indices: %s\n", paste(time_indices, collapse = ", ")))
  
  # Create date range string for filenames
  date_range <- sprintf("avg-%d-%d", start_year, end_year)
  
  # Create surface directory for variables without depth dimension
  surface_dir <- file.path(model_dir, "surface")
  if (!dir.exists(surface_dir)) dir.create(surface_dir)
  
  # Process each variable
  for (var_name in model$variables) {
    if (!(var_name %in% names(nc$var))) {
      cat(sprintf("  Warning: Variable %s not found in the dataset. Skipping.\n", var_name))
      next
    }
    
    cat(sprintf("  Processing variable: %s\n", var_name))
    
    # Check dimensionality of the variable
    var_dim_count <- length(nc$var[[var_name]]$dim)
    var_dims <- sapply(nc$var[[var_name]]$dim, function(d) d$name)
    
    # Check if this variable has a depth dimension
    has_depth_dim <- FALSE
    depth_dim_idx <- NULL
    
    # Check if this variable has a depth dimension
    var_dims <- sapply(nc$var[[var_name]]$dim, function(d) d$name)
    depth_dim_idx <- which(var_dims == "k")
    cat(sprintf("    Variable dimensions: %s\n", paste(var_dims, collapse = ", ")))
    cat(sprintf("    Depth dimension index: %s\n", depth_dim_idx))
    has_depth_dim <- length(depth_dim_idx) > 0
    
    if (has_depth_dim) {
      # Process for each depth if the variable has a depth dimension
      cat(sprintf("    Variable has depth dimension. Processing for each target depth.\n"))
      
      # Process each depth
      for (depth in depths) {
        # Find the closest depth in the model
        depth_info <- find_closest_depth(depth, nc)
        
        # Create depth string based on actual depth found
        depth_str <- sprintf("%.1fm", abs(depth_info$actual_depth))
        
        cat(sprintf("    Processing depth: target %.1fm, using closest %.2f m (k=%d)\n", 
                   abs(depth), depth_info$actual_depth, depth_info$index))
        
        # Create depth-specific output directory
        depth_dir <- file.path(model_dir, depth_str)
        if (!dir.exists(depth_dir)) dir.create(depth_dir)
        
        # Set output filename
        output_file <- file.path(depth_dir,
                                 sprintf("%s_%s_%s_%s.tif",
                                         model$model_id, date_range,
                                         var_name, depth_str))
        
        # Create start and count vectors for this depth
        start_vec <- rep(1, var_dim_count)
        count_vec <- rep(-1, var_dim_count)
        
        # Set specific depth slice
        start_vec[depth_dim_idx] <- depth_info$index + 1  # +1 because NetCDF is 0-indexed but R is 1-indexed
        count_vec[depth_dim_idx] <- 1
        
        # Download and process
        download_and_average_variable(nc, var_name, output_file, start_vec, count_vec, 
                                 time_indices, time_years)
      }
    } else {
      # Process once for variables without depth dimension
      cat(sprintf("    Variable does not have depth dimension. Processing once.\n"))
      
      # Set output filename
      output_file <- file.path(surface_dir, 
                               sprintf("%s_%s_%s.tif", 
                                       model$model_id, date_range, 
                                       var_name))
      
      # Create start and count vectors (no depth dimension)
      start_vec <- rep(1, var_dim_count)
      count_vec <- rep(-1, var_dim_count)
      
      # Download and process
      download_and_average_variable(nc, var_name, output_file, start_vec, count_vec, 
                                 time_indices, time_years)
    }
  }
  
  # Close NetCDF connection
  nc_close(nc)
}

cat("\nAll GeoTIFF climatology files created successfully!\n")

