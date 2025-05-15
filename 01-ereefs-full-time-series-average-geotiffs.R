# This script produces GeoTiff raster files corresponding to the average of
# the full time series of many of the common variable for the BGC and Hydro 
# model data for a set of depths. This creates a climatology estimate of 
# each variable allowing spatial patterns of average conditions to be identified. 
#
# This script downloads only a subset of the depth slices. A shallow depth that
# corresponds largely to just below the surface (3m) and a deeper depth that
# corresponds to the deeper coral areas (-17.75m).
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
# set of GeoTiff files saved in data/out/<Model><version>/<depth>/<ModelID>_<variable>_<depth>m.tif.
# Where <Model><version> are adjusted to replace " " with "_" and "." with "p".
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


library(ncdf4)   # For handling NetCDF files
library(raster)  # For creating and handling raster data
library(sf)      # For spatial operations

# Create output directory structure if it doesn't exist
output_dir <- "data/out"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Define the depths we want to extract (m)
depths <- c(-3, -18)

# Model definitions - contains all the information we need to access each model
models <- list(
#   list(
#     name = "GBR4 BGC baseline",
#     version = "3.1",
#     model_id = "GBR4_H2p0_B3p1_Cq3b_Dhnd",
#     url = "https://thredds.ereefs.aims.gov.au/thredds/dodsC/GBR4_H2p0_B3p1_Cq3b_Dhnd/annual.nc",
#     variables = c('TN', 'TP', 'DIN', 'DIP', 'Chl_a_sum', 'NO3', 'NH4', 'DOR_N', 
#                   'DOR_P', 'Secchi', 'PH', 'omega_ar', 'EFI'),
#     type = "bgc"
#   ),
#   list(
#     name = "GBR4 Hydro",
#     version = "2.0",
#     model_id = "gbr4_v2_hydro",
#     url = "https://thredds.ereefs.aims.gov.au/thredds/dodsC/gbr4_v2/annual.nc",
#     variables = c('eta', 'temp', 'salt', 'u', 'v', 'wspeed_u', 'wspeed_v'),
#     type = "hydro"
#   ),
  list(
    name = "GBR1 Hydro",
    version = "2.0",
    model_id = "gbr1_2.0_hydro",
    url = "https://thredds.ereefs.aims.gov.au/thredds/dodsC/gbr1_2.0/annual.nc",
    variables = c('eta', 'temp', 'salt', 'u', 'v', 'wspeed_u', 'wspeed_v'),
    type = "hydro"
  )
# Commented out for now as the model is not yet available
#   list(
#     name = "GBR4 Hydro",
#     version = "4.0",
#     model_id = "GBR4_H4p0_ABARRAr2_OBRAN2020_FG2Gv3_Dhnd",
#     url = "https://thredds.ereefs.aims.gov.au/thredds/dodsC/GBR4_H4p0_ABARRAr2_OBRAN2020_FG2Gv3_Dhnd/annual.nc",
#     variables = c('eta', 'temp', 'salt', 'u', 'v', 'wspeed_u', 'wspeed_v'),
#     type = "hydro"
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

# Function to sanitize model name for directory structure
sanitize_name <- function(name) {
  name <- gsub(" ", "_", name)
  name <- gsub("\\.", "p", name)
  return(name)
}

# Function to copy NetCDF attributes to the raster object
copy_nc_attributes_to_raster <- function(nc, var_name, r) {
  # Get all variable attributes
  attrs <- ncatt_get(nc, var_name)
  
  # Initialize the attributes list if it doesn't exist
  if (is.null(r@data@attributes) || length(r@data@attributes) == 0) {
    r@data@attributes <- list(list())
  }
  
  # Add these as metadata to the raster
  for (att_name in names(attrs)) {
    if (att_name != "_FillValue") {  # Skip fill value as it's handled differently in rasters
      r@data@attributes[[1]][[att_name]] <- attrs[[att_name]]
    }
  }
  
  # Add standard raster metadata if available
  if ("units" %in% names(attrs)) {
    r@data@unit <- attrs$units
  }
  if ("long_name" %in% names(attrs)) {
    r@data@names <- attrs$long_name
  } else {
    r@data@names <- var_name
  }
  
  return(r)
}

# Process each model
for (model in models) {
  model_dir_name <- paste0(sanitize_name(model$name), sanitize_name(model$version))
  cat(sprintf("\nProcessing model: %s %s\n", model$name, model$version))
  
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
  
  # Get coordinate information
  lon <- ncvar_get(nc, "longitude")
  lat <- ncvar_get(nc, "latitude")
  time <- ncvar_get(nc, "time")
  time_units <- ncatt_get(nc, "time", "units")$value
  
  # Process each depth
  for (depth in depths) {
    # First find the closest depth in the model using the zc variable
    depth_info <- find_closest_depth(depth, nc)
    
    # Now create depth_str based on the actual depth found
    depth_str <- sprintf("%.1fm", abs(depth_info$actual_depth))
    
    cat(sprintf("  Processing depth: target %.1fm, using closest %.2f m (k=%d)\n", 
                abs(depth), depth_info$actual_depth, depth_info$index))
    
    # Create depth-specific output directory
    depth_dir <- file.path(model_dir, depth_str)
    if (!dir.exists(depth_dir)) dir.create(depth_dir)
    
    # Process each variable
    for (var_name in model$variables) {
        
        # Get the number of time steps
        time_dim <- nc$dim$time$len
        cat(sprintf("    Processing variable: %s (%d time steps)\n", var_name, time_dim))
        
        # Check dimensionality of the variable
        var_dim_count <- length(nc$var[[var_name]]$dim)
        
        # Create appropriate filename based on whether variable has depth dimension
        if (var_dim_count == 4) {
            output_file <- file.path(depth_dir, 
                                sprintf("%s_%s_%s.tif", model$model_id, var_name, depth_str))
        } else {
            output_file <- file.path(depth_dir, 
                                sprintf("%s_%s.tif", model$model_id, var_name))
        }
        
        # Initialize arrays for accumulating values and counting non-NA cells
        first_slice <- TRUE
        sum_data <- NULL
        count_data <- NULL
        
        # Process one time step at a time. This is to limit the size of the OpenDAP
        # request, that can fail if it is too large.
        for (t in 1:time_dim) {
            cat(sprintf("\r    Processing time slice %d/%d...", t, time_dim))
            
            # Get data for this time slice
            if (var_dim_count == 4) {
                # 4D variable with depth dimension
                slice_data <- ncvar_get(
                    nc, var_name,
                    start = c(1, 1, depth_info$index + 1, t),
                    count = c(-1, -1, 1, 1)
                )
            } else if (var_dim_count == 3) {
                # 3D variable without depth dimension
                slice_data <- ncvar_get(
                    nc, var_name,
                    start = c(1, 1, t),
                    count = c(-1, -1, 1)
                )
            }
            
            # Remove singleton dimensions
            slice_data <- drop(slice_data)
            
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
        
        # Create raster
        r <- raster(
            t(mean_data),  # Transpose because raster expects rows=y, cols=x
            xmn = min(lon), xmx = max(lon),
            ymn = min(lat), ymx = max(lat),
            crs = "+proj=longlat +datum=WGS84"
        )
        
        # Copy metadata from NetCDF to raster
        r <- copy_nc_attributes_to_raster(nc, var_name, r)
        
        # Fix the vertical flip issue
        r <- flip(r, direction='y')
        
        # Save as GeoTIFF
        writeRaster(r, output_file, format = "GTiff", overwrite = TRUE)
        
        cat(sprintf("    Created GeoTIFF: %s\n", output_file))
    }
  }
  
  # Close NetCDF connection
  nc_close(nc)
}

cat("\nAll GeoTIFF climatology files created successfully!\n")

