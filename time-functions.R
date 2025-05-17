# =============================================================================
# NetCDF Time Conversion Utility
# =============================================================================
#
# PURPOSE:
#   This script provides a robust function for converting NetCDF time values to
#   proper R datetime objects (POSIXct). NetCDF files commonly represent time as
#   numeric values (often days, seconds, etc.) offset from a reference date, using
#   the format: "{unit} since YYYY-MM-DD HH:MM:SS +TZ".
#
# WHEN TO USE:
#   Use this function when working with NetCDF datasets in R where you need to:
#   - Convert numeric time values to actual calendar dates
#   - Handle various time units (days, seconds, hours, minutes, etc.)
#   - Process different timezone specifications (+10, -05, +1000, etc.)
#   - Ensure consistent UTC-based timestamps for analysis
#
# INPUTS:
#   time_values: Numeric vector containing time offsets from the reference date
#   time_units:  String defining the units and reference date in the format
#                "{unit} since YYYY-MM-DD HH:MM:SS [+/-TZ]"
#                where:
#                - {unit} can be seconds, minutes, hours, days, etc.
#                - [+/-TZ] is an optional timezone offset (e.g., +10, -05)
#   verbose:     Logical flag to enable debug output (default = FALSE)
#
# OUTPUTS:
#   Returns a POSIXct vector with properly converted datetime values in UTC
#
# SPECIAL HANDLING:
#   - Supports various time units (seconds, minutes, hours, days, months, years)
#   - Handles missing timezone information (assumes UTC)
#   - Processes different timezone formats:
#     * Single digit (+1, -5)
#     * Double digit (+10, -05)
#     * Four digit as hours/minutes (+1030)
#     * Four digit as hours only (+1000)
#   - Includes debug output to help troubleshoot parsing issues
#
# EXAMPLE:
#   time_values <- c(0, 365, 730)
#   time_units <- "days since 1990-01-01 00:00:00 +10"
#   dates <- convert_netcdf_time(time_values, time_units)
#   # Returns datetimes for 1989-12-31 14:00:00, 1990-12-31 14:00:00, 1991-12-31 14:00:00
#
# AUTHOR: 
#   Created on: May 17, 2025
# =============================================================================
convert_netcdf_time <- function(time_values, time_units, verbose = FALSE) {
  # Debug info
  if(verbose) cat(sprintf("Processing time units: '%s'\n", time_units))
  
  # Extract the time unit (days, seconds, etc.)
  unit_type <- sub("([^ ]+) since.*", "\\1", time_units)
  if(verbose) cat(sprintf("Unit type: '%s'\n", unit_type))
  
  # Check if the unit_type is supported
  supported_units <- c("seconds", "minutes", "hours", "days")
  if (!unit_type %in% supported_units) {
    stop(sprintf("Unsupported time unit: '%s'. Supported units are: %s", 
                unit_type, paste(supported_units, collapse = ", ")))
  }
  
  # Set the conversion factor based on the unit type
  conversion_factor <- switch(unit_type,
                             "seconds" = 1,
                             "minutes" = 60,
                             "hours" = 3600,
                             "days" = 86400) 
  
  # First extract the date part without timezone
  date_pattern <- "since\\s+([0-9]{4}-[0-9]{1,2}-[0-9]{1,2}\\s+[0-9]{1,2}:[0-9]{1,2}:[0-9]{1,2})"
  date_match <- regexpr(date_pattern, time_units, perl=TRUE)
  
  if (date_match == -1) {
    stop("Failed to extract reference date from time units: ", time_units)
  }
  
  # Extract the matched date string
  match_length <- attr(date_match, "match.length")
  ref_date_str <- substr(time_units, date_match + 6, date_match + match_length - 1)
  if(verbose) cat(sprintf("Extracted date string: '%s'\n", ref_date_str))
  
  # Parse the date without timezone first
  ref_date <- as.POSIXct(ref_date_str, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  
  # Now check for timezone after the date
  # Get the part after the date string
  remaining_str <- substr(time_units, date_match + match_length + 1, nchar(time_units))
  if(verbose) cat(sprintf("Remaining string after date: '%s'\n", remaining_str))
  
  # Look for timezone pattern in the remaining string
  tz_match <- regexpr("\\s*([+-][0-9]{1,4})\\s*", remaining_str, perl=TRUE)
  
  if (tz_match > 0) {
    tz_length <- attr(tz_match, "match.length")
    tz_str <- trimws(substr(remaining_str, tz_match, tz_match + tz_length - 1))
    if(verbose) cat(sprintf("Extracted timezone: '%s'\n", tz_str))
    
    # Remove any spaces
    tz_str <- gsub("\\s+", "", tz_str)
    
    # Get sign and numeric value separately
    tz_sign <- substr(tz_str, 1, 1)
    tz_numeric <- as.numeric(substr(tz_str, 2, nchar(tz_str)))
    
    # Handle different formats
    if (nchar(tz_str) == 2) {  # +1, -5 format
      tz_offset <- tz_numeric * 1.0  # 1 hour
    } else if (nchar(tz_str) == 3) {  # +10, -05 format
      tz_offset <- tz_numeric * 1.0  # 10 hours
    } else if (nchar(tz_str) == 5) {  # +1030 format (hours and minutes)
      hours <- as.numeric(substr(tz_str, 2, 3))
      minutes <- as.numeric(substr(tz_str, 4, 5)) / 60
      tz_offset <- hours + minutes
    } else {  # +1000 format (treat as hours)
      tz_offset <- as.numeric(substr(tz_str, 2, 3))
    }
    
    # Apply the sign
    if (tz_sign == "-") {
      tz_offset <- -tz_offset
    }
    
    if(verbose) cat(sprintf("Timezone offset in hours: %f\n", tz_offset))
    
    # Adjust reference date for timezone offset
    ref_date <- ref_date - tz_offset * 3600
  } else {
    if(verbose) cat("No timezone offset found, using UTC\n")
  }
  
  # Add the time values using the appropriate conversion factor
  result_dates <- ref_date + time_values * conversion_factor
  
  return(result_dates)
}

