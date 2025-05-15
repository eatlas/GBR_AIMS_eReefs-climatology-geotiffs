# Copyright Eric Lawrey, Australian Institute of Marine Science
# The script downloads the monthly BCG data for a set of variables for a specific
# depth, saving the downloaded data as a NetCDF file locally. This script also 
# copies over all the metadata from the source. 
#
# To limit the amount of data downloaded we are using the monthly average data and
# a single depth layer.
# 
# This script downloads eReefs BGC data from the AIMS THREDDS data service 
# using OpenDAP. It saves the data as a local NetCDF file for later processing.
# 
#
# One limitation of this script is that the copying of the metadata is slow because
# of the way the ncdf4 package works. It is not possible to copy all the metadata
# in one go. This is a known issue with the ncdf4 package.

library(ncdf4)  # Using only ncdf4 for all NetCDF operations

# We are using the OpenDAP service provided by AIMS to access the eReefs BGC data.
# The AIMS service provides access to monthly aggregated data which vastly reduces
# the amount of data we need to download. 
# Here we are using the GBR v3.1 base model, based on the Hydrographic model v2.0.
# In the near future (June 2025 onwards) we will need to update this to the 
# as yet unreleased GBR v4.x model.
url <- "https://thredds.ereefs.aims.gov.au/thredds/dodsC/GBR4_H2p0_B3p1_Cq3b_Dhnd/monthly.nc"

# We need to determine what is the depth index of the layer that most closely 
# matches our desired depth.
depth <- -3

# GBR4 depth-to-k lookup table
# Define the LUT as a numeric vector, where index = k + 1 (since R is 1-based)
# k = 0:16, so gbr4_depths[k+1] gives the depth at index k
gbr4_depths <- c(
  -145.0, -120.0, -103.0, -88.0, -73.0, -60.0, -49.0, -39.5,
  -31.0, -23.75, -17.75, -12.75, -8.8, -5.55, -3.0, -1.5, -0.5
)

# List of variables that we want to download from the BGC model. To determine the
# names of these variables or what they mean refer to the 'Output variables' section
# of  https://open-aims.github.io/ereefs-tutorials/tutorials/general/intro_to_ereefs.html
variables <- c('TN', 'TP', 'DIN', 'DIP', 'Chl_a_sum', 'NO3', 'NH4', 'DOR_N', 'DOR_P', 
               'PhyL_N', 'PhyL_NR', 'PhyS_N', 'PhyS_NR', 'Tricho_N', 'Tricho_NR')

# Filename of the output NetCDF file.
output_filename <- sprintf('GBR4_H2p0_B3p1_Cq3b_Dhnd_WQ_monthly_%sm.nc', abs(depth))

# Destination folder for the downloaded NetCDF file
destination_folder <- file.path("data","out","monthly")

# Function to copy the metadata associated with variables and dimensions. 
# This function is used to copy the metadata for the:
# - variables (e.g. TN, TP, etc.)
# - dimensions (e.g. time, latitude, longitude)
# - global attributes
copy_nc_attributes <- function(nc_source, nc_dest, src_var, dest_var, verbose = TRUE) {
  # The metadata is stored as a list of attributes in the NetCDF file.
  # These are key value pairs that describe the variable or dimension.
  all_atts <- ncatt_get(nc_source, src_var)
  
  # Copy each attribute
  for (att_name in names(all_atts)) {
    attval <- all_atts[[att_name]]
    # Note that this operation is slow. It takes about 5 sec per put
    # because it triggers a write to the NetCDF file.
    # This is a known issue with the ncdf4 package.
    ncatt_put(nc_dest, dest_var, att_name, attval)
    if (verbose) {
      cat(sprintf("Attribute %s: %s\n", att_name, attval))
    }
  }
}

if (!dir.exists(destination_folder)) dir.create(destination_folder, recursive = TRUE)

# Find the index of the closest depth
closest_idx <- which.min(abs(gbr4_depths - depth))
nearest_depth <- gbr4_depths[closest_idx]
k <- closest_idx - 1  # Convert R's 1-based index to k index

cat(sprintf("Requested depth: %.2f, using nearest LUT depth: %.2f (k=%d)\n", depth, nearest_depth, k))

cat("Opening the remote dataset...\n")
nc <- nc_open(url)

filename <- file.path(destination_folder, output_filename)

cat("Downloading data from OpenDAP service. This may take a while...\n")

# Print the details of the dimensions for one of the variables. This is useful
# debugging information because it shows the order of the dimensions and their sizes.
# The THREDDS OPeNDAP service (https://thredds.ereefs.aims.gov.au/thredds/dodsC/GBR4_H2p0_B3p1_Cq3b_Dhnd/daily.nc.html)
# lists the dimensions as:
# [time = 3073][k = 17][latitude = 723][longitude = 491]
# But the order when read using ncdf4 is:
# Variable 'TN' dimensions: [longitude = 491][latitude = 723][k = 17][time = 101]
varname <- variables[1]
dims <- nc$var[[varname]]$dim
dim_str <- paste0("[", 
    sapply(dims, function(d) sprintf("%s = %d", d$name, d$len)), 
    "]", collapse = "")
cat(sprintf("Variable '%s' dimensions: %s\n", varname, dim_str))

# Create a new NetCDF file with the same dimensions as the source
# Copy dimensions and coordinate variables
# Note: We drop the depth dimension (k) because we are only interested in a 
# single depth slice.
# Get coordinate values from source NetCDF first
lon_vals <- ncvar_get(nc, "longitude")
lat_vals <- ncvar_get(nc, "latitude")
time_vals <- ncvar_get(nc, "time")

lon_units <- ncatt_get(nc, "longitude", "units")$value
lat_units <- ncatt_get(nc, "latitude", "units")$value
time_units <- ncatt_get(nc, "time", "units")$value

# Define dimensions with their actual values and copied units
londim <- ncdim_def("longitude", lon_units, lon_vals)
latdim <- ncdim_def("latitude", lat_units, lat_vals)
timedim <- ncdim_def("time", time_units, time_vals)

# Create a list for data variables only (not coordinate variables)
varlist <- list()

# Define data variables
for (varname in variables) {
    # Create a variable definition.
    units <- "" # Blank for now, will be copied later by copy_nc_attributes
    dimensions <- list(londim, latdim, timedim)
    var_def <- ncvar_def(varname, units, dimensions, prec="float")
    varlist[[varname]] <- var_def
}

# Create the NetCDF file with data variables 
# (dimensions and coordinate variables are created automatically)
nc_new <- nc_create(filename, varlist)

# Copy coordinate metadata
for (coord in c("longitude", "latitude", "time")) {
    copy_nc_attributes(nc, nc_new, coord, coord)
}

# Download and add data variables
num_vars <- length(variables)
for (i in seq_along(variables)) {
    var <- variables[i]
    cat(sprintf("(%d/%d) Fetching data for %s at depth=%s... (This will take ~1 min)\n", 
                i, num_vars, var, depth))
    # Download only the required depth slice
    data_k <- ncvar_get(
        nc, var,
        start = c(1, 1, k + 1, 1),
        count = c(length(lon_vals), length(lat_vals), 1, length(time_vals))
    )

    # Remove the singleton k dimension and reshape for [lon, lat, time]
    data_k <- array(data_k, dim = c(length(lon_vals), length(lat_vals), length(time_vals)))

    # Save the data
    ncvar_put(nc_new, var, data_k)
    
    # Copy the metadata
    copy_nc_attributes(nc, nc_new, var, var)
    
    cat(sprintf("Data for %s saved successfully!\n", var))
}

# Copy global attributes
copy_nc_attributes(nc, nc_new, 0, "NC_GLOBAL")

nc_close(nc_new)
nc_close(nc)
cat("All data extraction and saving completed!\n")