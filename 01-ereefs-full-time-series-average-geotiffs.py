"""
This script produces GeoTiff raster files corresponding to the average of
the full time series of many common variables for the eReefs BGC and Hydro 
model data across multiple depths. This creates a climatology estimate of 
each variable allowing spatial patterns of average conditions to be identified.

The script processes four standard depths: 3m, 9m, 18m, and 39m below the surface.
Each model is processed using predefined date ranges to ensure only complete years 
are included. Variables without a depth dimension are processed once and saved 
in a 'surface' directory.

Data is obtained from the AIMS THREDDS server via OPeNDAP, aggregated across
all full years, and saved as GeoTiff files with appropriate metadata.
"""

import os
import numpy as np
import xarray as xr
import rasterio
from rasterio.transform import from_origin
from pathlib import Path
import pandas as pd  # For timestamp handling
import warnings
warnings.filterwarnings('ignore')  # Suppress warnings from xarray/netcdf4

# Create output directory structure
output_dir = Path("data/out-py")
output_dir.mkdir(parents=True, exist_ok=True)

# Define depths to extract (m)
depths = [-3, -9, -18, -39]

# Model definitions
models = [
    {
        "name": "GBR4-BGC3p1-base",
        "model_id": "GBR4_H2p0_B3p1_Cq3b_Dhnd",
        "url": "https://thredds.ereefs.aims.gov.au/thredds/dodsC/GBR4_H2p0_B3p1_Cq3b_Dhnd/annual.nc",
        "variables": ['TN', 'TP', 'DIN', 'DIP', 'Chl_a_sum', 'NO3', 'NH4', 'DOR_N', 
                    'DOR_P', 'Secchi', 'PH', 'omega_ar', 'EFI'],
        "type": "bgc",
        "start_year": 2011,
        "end_year": 2018
    },
    {
        "name": "GBR4-H2p0",
        "model_id": "gbr4_v2_hydro",
        "url": "https://thredds.ereefs.aims.gov.au/thredds/dodsC/gbr4_v2/annual.nc",
        "variables": ['eta', 'temp', 'salt', 'u', 'v', 'wspeed_u', 'wspeed_v'],
        "type": "hydro",
        "start_year": 2011,
        "end_year": 2023
    },
    {
        "name": "GBR1-H2p0",
        "model_id": "gbr1_2.0_hydro",
        "url": "https://thredds.ereefs.aims.gov.au/thredds/dodsC/gbr1_2.0/annual.nc",
        "variables": ['eta', 'temp', 'salt', 'u', 'v', 'wspeed_u', 'wspeed_v'],
        "type": "hydro",
        "start_year": 2015,
        "end_year": 2023
    }
]

def find_closest_depth(target_depth, dataset):
    """Find the closest depth index in the model to our target depth."""
    if 'zc' in dataset:
        depth_values = dataset['zc'].values
        idx = np.argmin(np.abs(depth_values - target_depth))
        return {
            "index": idx,
            "actual_depth": float(depth_values[idx])
        }
    else:
        raise ValueError("'zc' variable not found in the dataset. Cannot determine depth levels.")

def create_geotiff(data_array, output_file):
    """Create a GeoTIFF file from an xarray DataArray."""
    # Get coordinate information
    lon = data_array.longitude.values
    lat = data_array.latitude.values
    
    # Calculate pixel size and origin
    pixel_width = (lon[-1] - lon[0]) / (len(lon) - 1)
    pixel_height = (lat[-1] - lat[0]) / (len(lat) - 1)
    
    # eReefs data is latitudes ascending (south to north)

    # In rasterio, origin needs to be the upper-left (northwest) corner of the upper-left pixel
    # Need to adjust by half pixel to account for center vs corner registration

    # Latitude values go from south to north, so lat[-1] is the northernmost point
    # Adjust origin by half a pixel width/height to account for center vs corner registration
    transform = from_origin(
        lon[0] - (pixel_width/2),  # Shift half pixel west
        lat[-1] + (abs(pixel_height)/2),  # Shift half pixel north
        pixel_width, 
        abs(pixel_height)
    )
    # Need to flip the data because GeoTIFF expects north-up orientation
    data_values = np.flipud(data_array.values)
    
    
    # Remove any singleton dimensions
    data_values = np.squeeze(data_values)
    
    # Create GeoTIFF
    with rasterio.open(
        output_file,
        'w',
        driver='GTiff',
        height=data_values.shape[0],
        width=data_values.shape[1],
        count=1,
        dtype=data_values.dtype,
        crs='+proj=longlat +datum=WGS84',
        transform=transform,
        nodata=np.nan
    ) as dst:
        dst.write(data_values, 1)
        
        # Add metadata (variable attributes)
        for attr_name, attr_value in data_array.attrs.items():
            if isinstance(attr_value, (str, int, float)):
                dst.update_tags(**{attr_name: str(attr_value)})
    
    print(f"    Created GeoTIFF: {output_file}")

def download_and_average_variable(dataset, var_name, output_file, start_index, end_index, k_index=None):
    """Download, average, and save a variable as GeoTIFF using time indices.
    
    Processes one time slice at a time to conserve memory and provide detailed statistics.
    """
    time_indices = np.arange(start_index, end_index + 1)
    
    # Get QLD times for reporting
    qld_time_values = convert_to_qld_time(dataset.time.values)
    qld_years = [pd.Timestamp(t).year for t in qld_time_values]
    
    print(f"\n    Time slice averages for {var_name}:")
    
    # Initialize a list to store individual time slices for later averaging
    time_slices = []
    
    # Process each time slice individually
    for idx in time_indices:
        # Access just this time slice
        if k_index is not None:
            # For 3D variables (with depth)
            time_slice = dataset[var_name].isel(time=idx, k=k_index)
        else:
            # For 2D variables (no depth)
            time_slice = dataset[var_name].isel(time=idx)
        
        # Calculate average (mean) of this time slice, ignoring NaN values
        slice_avg = float(time_slice.mean(skipna=True).values)
        print(f"      Time index {idx} (year {qld_years[idx]}): {slice_avg:.6f}")
        
        # Store the time slice for later averaging
        time_slices.append(time_slice)
    
    # Combine all time slices along a new time dimension
    # This step temporarily uses more memory but ensures proper weighted averaging
    combined_data = xr.concat(time_slices, dim='time')
    
    # Calculate mean across time dimension
    mean_data = combined_data.mean(dim='time', skipna=True)
    
    # Print the overall mean for validation
    overall_mean = float(mean_data.mean(skipna=True).values)
    print(f"\n    Overall mean: {overall_mean:.6f}")
    
    # Create GeoTIFF
    create_geotiff(mean_data, output_file)
    
    return mean_data


def convert_to_qld_time(utc_times):
    """Convert UTC time values to Queensland local time (UTC+10)"""
    return np.array([pd.Timestamp(t) + pd.Timedelta(hours=10) for t in utc_times])

# Process each model
for model in models:
    model_dir_name = model["name"]
    print(f"\nProcessing model: {model['name']}")
    
    # Create model output directory
    model_dir = output_dir / model_dir_name
    model_dir.mkdir(exist_ok=True)
    
    # Create surface directory for variables without depth dimension
    surface_dir = model_dir / "surface"
    surface_dir.mkdir(exist_ok=True)
    
    
    try:
        # Open dataset with chunking to optimize memory usage
        ds = xr.open_dataset(model["url"], chunks={'time': 1})
        
        # Get years using QLD time for naming
        qld_time_values = convert_to_qld_time(ds.time.values)
        time_years = [pd.Timestamp(t).year for t in qld_time_values]

        # Determine the start and end for the start and end years based on time_years
        start_year = model["start_year"]
        end_year = model["end_year"]
        start_index = time_years.index(start_year) if start_year in time_years else -1
        end_index = time_years.index(end_year) if end_year in time_years else -1
        if start_index == -1 or end_index == -1:
            msg = f"Start year {start_year} or end year " \
                f"{end_year} not found in the dataset. Available years: {time_years}"
            raise ValueError(msg) 
        
        
        # Create date range string for filenames
        date_range = f"{start_year}-{end_year}"
        
        print(f"  Using time indices: {start_index} to {end_index} for years {start_year}-{end_year}")
        
        # Process each variable
        for var_name in model["variables"]:
            if var_name not in ds:
                print(f"  Warning: Variable {var_name} not found in the dataset. Skipping.")
                continue
            
            print(f"  Processing variable: {var_name}")
            
            # Check if variable has depth dimension
            has_depth_dim = 'k' in ds[var_name].dims
            print(f"    Variable dimensions: {', '.join(ds[var_name].dims)}")
            print(f"    Has depth dimension: {has_depth_dim}")
            
            if has_depth_dim:
                # Process for each depth
                print(f"    Variable has depth dimension. Processing for each target depth.")
                
                for depth in depths:
                    # Find closest depth in model
                    depth_info = find_closest_depth(depth, ds)
                    depth_str = f"{abs(depth_info['actual_depth']):.1f}m"
                    
                    print(f"    Processing depth: target {abs(depth):.1f}m, using closest "
                          f"{depth_info['actual_depth']:.2f}m (k={depth_info['index']})")
                    
                    # Create depth-specific output directory
                    depth_dir = model_dir / depth_str
                    depth_dir.mkdir(exist_ok=True)
                    
                    # Set output filename
                    output_file = depth_dir / f"{model['model_id']}_{date_range}_{var_name}_{depth_str}.tif"
                    
                    # Download and process
                    download_and_average_variable(
                        ds, var_name, output_file, start_index, end_index, k_index=depth_info['index']
                    )
            else:
                # Process once for variables without depth
                print(f"    Variable does not have depth dimension. Processing once.")
                
                # Set output filename
                output_file = surface_dir / f"{model['model_id']}_{date_range}_{var_name}.tif"
                
                # Download and process
                download_and_average_variable(
                    ds, var_name, output_file, start_index, end_index
                )
                
        # Close dataset
        ds.close()
        
    except Exception as e:
        print(f"  Error processing model {model['name']}: {str(e)}")
        exit(1)

print("\nAll GeoTIFF climatology files created successfully!")