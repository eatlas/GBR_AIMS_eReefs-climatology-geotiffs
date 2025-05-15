# eReefs data services - raw model data verse derived regridded aggregate data
eReefs data is provided in both its original raw modeled format and in regridded form with temporal aggregation:

## Raw model data: 
This data has the highest temporal resolution: hourly for Hydrodynamic models and daily for BioGeoChemical models. The raw model data also provides the full range of simulated depths from shallow -0.5 m through to the abyss -3880 m. The two complications with using the raw model data is that the data uses curvilinear grids, which makes it more difficult to work with, and the data files can be very large if the processing requires a long time series.
The raw model data is provided on the NCI Thredds data service from https://thredds.nci.org.au/thredds/catalog/catalogs/fx3/catalog.html The raw model data is developed by CSIRO then published on NCI.

## Regridded aggregate data:  
The model data is also available in a regridded form using a regular rectangular grid, which is compatible with most GIS tools and programming libraries. This collection is derived from the raw model data and is available from the AIMS eReefs Thredds data service. This collection also includes temporal aggregations of the data, providing daily averages from hourly data, and monthly and annual averages. The daily hydrodynamic data is an average of the raw hourly data, while the daily BGC data is a time snapshot, since there is no raw hourly BGC data. If the research question can work with averages over time then using the temporal aggregate data is much easier as the data size is much smaller. For a particular variable and depth the monthly aggregate is ~30x  smaller than the raw BGC data and ~30x24x smaller hydrodynamic data. There is no regridded version of the hourly hydrodynamic data due to the high storage costs associated such a dataset. Additionally the regridded dataset also only contains the top half of the modelled depth layers. This was to save on store costs. If you need depths greater than 145 m then you will need to use the raw model data. If you need to have hourly hydrodynamic data then you will need to use the raw model data. 
The regridded aggregate data is derived from the raw model data. This is processed by AIMS and published to a Thredds data service maintained by AIMS. These data products are available from: https://thredds.ereefs.aims.gov.au/thredds/catalog/catalog.html


The BioGeoChemical model is available in several scenarios. Each corresponding to different land run off conditions. The Baseline scenario corresponds to  the best estimate of historical conditions. The pre-industrial scenario corresponds to estimate river loads if the land use matched pre-industrial conditions. In this scenario the weather conditions match the baseline scenario. This means that changes seen in the marine environment between the two scenarios is because of the sediment and nutrient loads and not different weather conditions. The reduced scenario corresponds to river loads that correspond to water quality improvement targets.
 
In this script we want to extract estimates for historical conditions on the GBR and thus we want to use the baseline version of the BGC model.

The AIMS Thredds service provides multiple ways of accessing the data. The data is originally organised into a time series of NetCDF files. These NetCDF files are available as individual files, available via the HTTPServer service from THREDDS, or as a virtual dataset that stitches all the files into a single data source that represents the entire time series. This is then available as an OPeNDAP service. 

If you want to download a complete copy of the eReefs model data (all variables, all depths) then downloading from the HTTPServer service is the most efficient  approach. The HTTPServer service can be used to download individual files using a web browser.
The OPeNDAP service is intended to be access programmatically using libraries that understand the OPeNDAP protocol. OPeNDAP allows slices of data to be accessed, often eliminating the need to download all the data. You can just download  variables, time slices and depth slices of interest. Libraries in R and Python make performing these operations easy.
Each original NetCDF model file has both an OPeNDAP and HTTPServer service. It is thus possible to connect your code to use OPeNDAP to access a subset of an individual file. This is fine if you are looking at a single time point, but if you are looking at a longer time series, connecting to the data end point for each files is fiddly. For this is it much better to use the virtual dataset. It allows the entire time series to be accessed from a single OPeNDAP service end point.

In this script we want to create a climatology of the entire time series of data for a set of variables for a single depth slice. Since we only want a subset of variables and depths it is much easier using the virtual dataset OPeNDAP service.

# Depth layers
eReefs has horizontal depth layers (Z-Levels) that correspond to fixed depths regardless of the underlying water depth. The spacing between the depth layers are not equally spaced. There are more layers representing the shallow portions than deep water. Each depth slice is approximately 20 - 50% deeper than the next shallowest layer. The depth is represented by an index k into the depth dimension (zc), where a k of 0 corresponds to the deepest depth slice. For older versions of eReefs, version 2.0 of the hydro model and v3.1 of the BGC model, there is a slight difference in the depths of each slices between the GBR1 model grid and the GBR4 grid. 
In more recent versions 4.x the depth slices have been standardised across the models. This allows the GBR1 model to be more easily nested into the GBR4 model grid. The list of depths associated with a model can be obtains by looking at the metadata available in the OPeNDAP Dataset Access Form. For the BGC GBR4 v3.1 model from the regridded aggregate AIMS Thredds service this corresponds to https://thredds.ereefs.aims.gov.au/thredds/dodsC/GBR4_H2p0_B3p1_Cq3b_Dhnd/daily.nc.html. On this page list of depths can be obtained by checking the zc variable then clicking on the 'Get ASCII' Action button. This will return some short metadata and the list of depths:
```
Dataset {
  Float64 zc[k = 17];
} GBR4_H2p0_B3p1_Cq3b_Dhnd/daily.nc;
---------------------------------------------
  zc[17]
-145.0, -120.0, -103.0, -88.0, -73.0, -60.0, -49.0, -39.5, -31.0, -23.75, -17.75, -12.75, -8.8, -5.55, -3.0, -1.5, -0.5
```
Here k of 0 corresponds to -145 m depth and k of 16 corresponds to -0.5 m.