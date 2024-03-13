import os
import s3fs
import boto3
import fsspec
import pandas as pd
import numpy as np
import xarray as xr
import zarr
import glob
import rasterio
import pyproj
import geopandas
import matplotlib.pyplot as plt

# # set up Dask 
# import dask
# from dask.distributed import Client
# from dask.distributed import progress
# try:
#     print(client.dashboard_link)
# except:    
#     client = Client(n_workers=24, threads_per_worker=1, memory_limit='2GB') 


# functions
def get_nwm_id(geo_lat, geo_lon, ds):
    '''
    This function determines the NWM ID linked to a nearby reach based on a specified geographic location.
    input: 
        geo_lat: latitude of the geo-location
        geo_lon: longitude of the geo-location
        ds: the NWM dataset that includes streamflow data
    output: 
        NWM reach/feature ID
    '''
    # Calculate the Euclidean distance to find the closest point
    distances = np.sqrt((ds.latitude - geo_lat)**2 + (ds.longitude - geo_lon)**2)
    
    # Find the index of the minimum distance
    min_distance_index = np.unravel_index(distances.argmin(), distances.shape)

    # Extract the subset based on the indices
    subset = ds.isel(feature_id=min_distance_index[0].item()).compute()
    
    return subset.feature_id.values
    
    
    
def get_nwm_q(nwm_reach_id, start_date, end_date, output_path):
    '''
    This function subsets the NWM streamflow data for a specified stream reach and period.
    input:
        nwm_reach_id: The unique ID of stream reaches as defined in the NWM.
        start_date: The start date of the period of interest, formatted as 'yyyy-mm-dd'
        end_date: The end date of the period of interest, formatted as 'yyyy-mm-dd'
        
    output:
        output_path: The complete path to a CSV file to save subset results.
    '''
    print('Loading the metadata')
    conus_bucket_url = 's3://noaa-nwm-retrospective-3-0-pds/CONUS/zarr/chrtout.zarr'
    ds = xr.open_zarr(fsspec.get_mapper(conus_bucket_url, anon=True), consolidated=True)
    
    print('Subset started ...')
    start_date = pd.to_datetime(start_date)
    end_date = pd.to_datetime(end_date)
    subset = ds.sel(feature_id=int(nwm_reach_id)).sel(time=slice(start_date, end_date)).persist()
    subset.compute()
    subset.close()
    print('Subset finished!')
    
    print('Saving results in a csv file ...')
    subset_df = subset.streamflow.to_dataframe()
    subset_df.to_csv(output_path)
    print('Subset is completed.')
    
    return subset_df