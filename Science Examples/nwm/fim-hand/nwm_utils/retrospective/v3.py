#!/usr/bin/env python3

"""
Purpose: Utility functions for working with NWM retrospective data
Authors: Tony Castronova <acastronova@cuahsi.org>
         Irene Garousi-Nejad <igarousi@cuahsi.org>
Last Modified: March 26, 2024
"""

import numpy as np
import xarray as xr
import pandas as pd
from scipy import interpolate
import fsspec



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

    print('Loading the metadata...', end='')
    conus_bucket_url = 's3://noaa-nwm-retrospective-3-0-pds/CONUS/zarr/chrtout.zarr'
    ds = xr.open_zarr(fsspec.get_mapper(conus_bucket_url, anon=True), consolidated=True)
    print('done')
    
    print('Subset started...', end='')
    start_date = pd.to_datetime(start_date)
    end_date = pd.to_datetime(end_date)
    subset = ds.sel(feature_id=int(nwm_reach_id)).sel(time=slice(start_date, end_date))
    print('done.')
    
    print('Saving results in a csv file...', end='')
    subset_df = subset.streamflow.to_dataframe()
    subset_df.rename({'streamflow': f'streamflow_{ds.streamflow.attrs["units"]}'}, axis='columns', inplace=True)
    subset_df.to_csv(output_path)
    print('done')
    
    return subset_df

def interpolate_y(df: pd.DataFrame,
                  x_column: str,
                  y_column: str,
                  x_value: float) -> float:
    """
    Performs 1D interpolation on two columns of a Dataframe.

    Parameters
    ==========
    df: pandas.DataFrame
        DataFrame containing data that will be used in the interpolation.
    x_column: str
        Name of the column that represents the X-axis data.
    y_column: str
        Name of the column that represents the Y-axis data.
    x_value: float
        Numeric X-axis value for which to interpolate Y-axis data. Returns 
        -9999 if the interpolation fails to resolve.

    Returns
    =======
    y_value: float
        Numeric Y-axis value corresponding to the input X-axis value.

    """
    # Sort the DataFrame by the 'x' column to ensure interpolation works correctly
    df_sorted = df.sort_values(by=x_column)
    
    # Check if the x_value is within the range of the DataFrame
    if x_value < df_sorted[x_column].min() or x_value > df_sorted[x_column].max():
        return -9999  # x_value is out of range, cannot interpolate
    
    # Perform linear interpolation
    f = interpolate.interp1d(df_sorted[x_column], df_sorted[y_column], kind='linear')
    
    # Return the interpolated y value for the given x value
    return f(x_value)

