library(ncdf4)   # For handling NetCDF files
nc <- nc_open("https://thredds.ereefs.aims.gov.au/thredds/dodsC/GBR4_H2p0_B3p1_Cq3b_Dhnd/annual.nc")
time_var <- ncvar_get(nc, "time")
cat(sprintf("Time variable values: %s\n", paste(time_var, collapse = ", ")))


time_units <- ncatt_get(nc, "time", "units")$value
cat(sprintf("Time units: %s\n", time_units))

cat("---- Robust time conversion ----\n")
# This is an example of how to robustly handle the time conversion
# For this we must consider the time units and the timezone
# Most of the complexity is in the convert_netcdf_time function.
source("time-functions.R", local = TRUE)
# Convert to timezone-aware date-time objects based on the time units
date_time <- convert_netcdf_time(time_var, time_units)
cat(sprintf("Converted date-time values: %s\n", paste(date_time, collapse = ", ")))
# Convert to local time (UTC+10)
local_times <- format(date_time, tz="Australia/Brisbane")
# Extract years as numeric values
years <- as.numeric(format(as.POSIXct(local_times, tz="Australia/Brisbane"), "%Y"))
cat("Local times (UTC+10):", paste(local_times, collapse=", "), "\n")
cat("Years:", paste(years, collapse=", "), "\n")

cat("---- Simplified time conversion ----\n")
# If we want a simpler, slightly less robust version we can work off
# the fact that we know the units are days and the data is in UTC+10
# and the start date is 1990-01-01. We bake in these assumptions
# and simplify have a test to make sure these assumptions are correct.
assumed_time_units <- "days since 1990-01-01 00:00:00 +10"

# Fail if time units doesn't match
if (time_units != assumed_time_units) {
    stop(sprintf("Assumed time units '%s' do not match actual time units '%s'", 
                 assumed_time_units, time_units))
}

# Since we only want the year value and we know the units are days 
# and the day values are in local time we can just adjust the starting
# date by the time_var (in days).
year = round(1990 + time_var/365.25)
cat(sprintf("Years: %s\n", paste(year, collapse = ", ")))


