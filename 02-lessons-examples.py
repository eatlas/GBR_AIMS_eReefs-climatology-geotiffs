import xarray as xr
import pandas as pd

ds = xr.open_dataset("https://thredds.ereefs.aims.gov.au/thredds/dodsC/GBR4_H2p0_B3p1_Cq3b_Dhnd/annual.nc")
print(f"Dataset time values (UTC): {ds.time.values}")

# Convert to local time
local_time = [pd.Timestamp(t) + pd.Timedelta(hours=10) for t in ds.time.values]

# Extract the years
years = [pd.Timestamp(t).year for t in local_time]
print(f"Years: {years}")
