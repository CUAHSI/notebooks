#!/usr/bin/env python3

"""
Purpose: Utility functions for working with NOAA OWP FIM data
Authors: Tony Castronova <acastronova@cuahsi.org>
         Irene Garousi-Nejad <igarousi@cuahsi.org>
Last Modified: March 12, 2024
"""


import numpy as np
import pandas as pd
from typing import Dict
from pathlib import Path
from scipy import interpolate
import xarray as xr
import fsspec
from osgeo import gdal
import rasterio
import numpy.ma as ma

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
    subset_df.rename({'streamflow': f'streamflow_{ds.streamflow.attrs["units"]}'}, axis='columns', inplace=True)
    subset_df.to_csv(output_path)
    print('Subset is completed.')
    
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

def compute_stage(df: pd.DataFrame,
                  hydro_id: int,
                  flow_cms: float) -> Dict[int, float]:
    """
    Computes river stage from a rating curve given streamflow in cms.

    Parameters
    ==========
    df: pandas.DataFrame
        DataFrame containing the stage and discharge values of the rating curve.
        This must contain the following columns: HydroID, stage, discharge_cms.
    hydro_id: int
        Identifier for the reach for which to compute stage.
    flow_cms: float
        Streamflow to convert into river stage.

    Returns
    =======
    Dict [int, float]
        A dictionary containing computed stage and its associated hydroid

    """

    
    # look up rating curve for this hydroid
    rating_curve = df.loc[df.HydroID == hydro_id, ['stage', 'discharge_cms']]

    # interpolate using the provided flow rate
    interpolated_stage = interpolate_y(rating_curve, 'discharge_cms', 'stage', flow_cms) 

    return {hydro_id: float(interpolated_stage)}

def get_stage_for_nhd_reach(nhd_feature_id: int,
                            flow_cms: float,
                            hydrotable: Path = Path('./hydroTable_0.csv')) -> Dict[int, float]:
    """
    Retrieves stage for a given NHD reach and input streamflow. The stage for
    the most downstream HydroID is computed using FIM rating curves.

    Parameters
    ==========
    nhd_feature_id: int
        NHD feature identifier.
    flow_cms: float
        Streamflow for the reach in cubic meters per second
    hydrotable: pathlib.Path
        Path to the FIM hydrotable.csv file containing rating curve data.

    Returns
    =======
    Dict [int, float]
        Dictionary containing a single NWM hydro identifier and 
        its corresponding stage.

    """
    
    # load hydrotable_0
    # we don't need all of the columns in this csv
    hydro_df = pd.read_csv(hydrotable,
                           usecols=['HydroID', 'NextDownID', 'feature_id',
                                    'stage', 'discharge_cms'])
    
    # select features that match nhd_feature_id
    d = hydro_df.loc[hydro_df.feature_id==nhd_feature_id]
    
    # get unique combos of HydroID and NextDownID 
    from_nodes, to_nodes = np.unique(d[['HydroID', 'NextDownID']], axis=0).T
    
    # find the most downstream reach
    # computed as the 'From_Node' corresponding to the 'To_Node' that
    # that does not exist in the list of 'From_Node'
    downstream_hydroid = None
    for i, to_node in enumerate(to_nodes):
        if to_node not in from_nodes:
            # this is a node beyond our NHD reach
            # save the corresponding "from_node" id
            downstream_hydroid = from_nodes[i]
            
    # TODO: insert error handling here
    if downstream_hydroid == None:
        print('Something went wrong :( ')
        return {-9999: -9999.}
    
    interpolated_stage = compute_stage(d, downstream_hydroid, flow_cms)
    
    return interpolated_stage


def get_stage_for_all_hydroids_in_reach(nhd_feature_id: int,
                                        flow_cms: float,
                                        hydrotable: Path = Path('./hydroTable_0.csv')) -> Dict[int, float]:
    """
    Retrieves stage for all NWM HydroIDs given an NHD reach and input streamflow.
    The stage is computed using FIM rating curves.

    Parameters
    ==========
    nhd_feature_id: int
        NHD feature identifier.
    flow_cms: float
        Streamflow for the reach in cubic meters per second
    hydrotable: pathlib.Path
        Path to the FIM hydrotable.csv file containing rating curve data.

    Returns
    =======
    Dict [int, float]
        Dictionary containing one or more NWM hydro identifiers and 
        their corresponding stages.

    """
    
    # load hydrotable_0
    # we don't need all of the columns in this csv
    hydro_df = pd.read_csv(hydrotable,
                           usecols=['HydroID', 'NextDownID', 'feature_id',
                                    'stage', 'discharge_cms'])
    
    # select features that match nhd_feature_id
    d = hydro_df.loc[hydro_df.feature_id==nhd_feature_id]
    
    # get unique combos of HydroID and NextDownID 
    hydro_ids = np.unique(d.HydroID)
    
    interpolated_stages = {}
    for hydro_id in hydro_ids:  
        interpolated_stages.update(compute_stage(d, hydro_id, flow_cms))
        
    # return interpolated stage
    return interpolated_stages

def create_masked_stage_raster(input_raster, cumulative_mask_array, hydroID, stage):
    # Open the input raster file
    ds = gdal.Open(input_raster)
    if ds is None:
        print("Could not open input raster file")
        return
    
    # Read raster band as a numpy array
    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray()
    
    # Create a mask for the target value
    mask = np.where(arr == hydroID, stage, -99999)
    
    # Update the cumulative mask array
    if cumulative_mask_array is None:
        cumulative_mask_array = mask
    else:
        cumulative_mask_array = np.where(mask != -99999, mask, cumulative_mask_array)
    
    # Close the input dataset
    ds = None
    
    return cumulative_mask_array

def create_masked_stage_rasters(input_raster, output_masked_raster, res):
    # Initialize an empty cumulative mask array
    cumulative_mask_array = None
    
    # Open the input raster file to get its properties
    input_ds = gdal.Open(input_raster)
    if input_ds is None:
        print("Could not open input raster file")
        return
    
    # Loop over the dictionary 
    for hydroID, stage in res.items():
        print(f"HydroID: {hydroID}, Stage (m): {stage}")
        cumulative_mask_array = create_masked_stage_raster(input_raster, cumulative_mask_array, hydroID, stage)
    
    # Create the masked raster using the final cumulative mask array
    if cumulative_mask_array is not None:
        # Create a new raster with the cumulative mask
        driver = gdal.GetDriverByName("GTiff")
        out_ds = driver.Create(output_masked_raster, input_ds.RasterXSize, input_ds.RasterYSize, 1, gdal.GDT_Float32)
        out_ds.SetProjection(input_ds.GetProjection())
        out_ds.SetGeoTransform(input_ds.GetGeoTransform())
        out_band = out_ds.GetRasterBand(1)
        out_band.WriteArray(cumulative_mask_array)
        
        # Close the output dataset
        out_ds = None
        
        print("Masked raster created successfully.")
    else:
        print("No valid mask data to write.")
        

def create_masked_flood_inundation(hand_raster, stage_raster, output_raster):
    with rasterio.open(hand_raster) as src1, rasterio.open(stage_raster) as src2:
        hand = src1.read(1) #, masked=True
        hand = hand.astype(float)
        hand = ma.masked_equal(hand, -999999.0)
        print(hand.min(), hand.max())
        stage = src2.read(1) #, masked=True
        stage = stage.astype(float)
        stage = ma.masked_equal(stage, -99999.0)
        print(stage.min(), stage.max())
        transform = src2.transform

        # Perform subtraction and handle NaN values 
        hand_temp = src1.read(1, masked=True)
        stage_temp = src2.read(1, masked=True)
        result = np.where(stage_temp.mask, np.nan, hand_temp - stage_temp)
        result = result.astype(float)
        mask_condition = (result>0) & (result <-500) # (result < -500) | (result > 500) returns all the area  # (result < 0) | (result > 500) returns non flooded area within the domain
        masked_result = np.ma.masked_array(result, mask=mask_condition)
        # masked_result = ma.masked_equal(result, -999999.0) # if turned on -> min and max are not printed out correctly
        print(masked_result.min(), masked_result.max())
        # threshold = 0.00001  
        # result = np.where(result < threshold, 1, 0)
        
        # Set x_min and x_max values in the transform 
        masked_result = ma.masked_equal(result, -999999.0) # required to make sure the transform is built correctly.
        non_nan_indices = np.where(~masked_result.mask)
        xmin, ymin = transform * (non_nan_indices[1].min(), non_nan_indices[0].min())
        xmax, ymax = transform * (non_nan_indices[1].max(), non_nan_indices[0].max())
        print(xmin, xmax)
        transform = rasterio.transform.from_bounds(xmin, src2.bounds.bottom, xmax, src2.bounds.top, src2.width, src2.height)

        # Define the output raster metadata based on one of the input rasters
        profile = src1.profile  # Using profile from the hand raster
        # profile.update(dtype=rasterio.float32, nodata=np.nan)  # Update data type and nodata value
        profile['transform'] = transform  # add specified xmin and max to results
        
    # Write the result to a new raster file
    with rasterio.open(output_raster, 'w', **profile) as dst:
        masked_result = np.ma.masked_array(result, mask=mask_condition)
        dst.write(masked_result.astype(rasterio.float32), 1)
        
    src1.close()
    src2.close()
    

if __name__ == "__main__":
    nhd_feature_id=7897865
    cms = 2000

    print("\nTest NHD+ reach stage:\n")
    res = get_stage_for_nhd_reach(nhd_feature_id, cms)
    print(f'{"NHD FeatureID":^15}  {"HydroID":^15}  {"Input CMS":^15}  {"Output Stage":^15}')
    print(f'{"===============":15}  {"===============":15}  {"===============":15}  {"===============":15}')
    for k,v in res.items():
        print(f'{nhd_feature_id:<15}  {k:<15}  {CMS:<15}  {v:<15}')
    
    print("\nTest HydroID reach stages:\n")
    res = get_stage_for_all_hydroids_in_reach(nhd_feature_id, CMS)
    print(f'{"NHD FeatureID":^15}  {"HydroID":^15}  {"Input CMS":^15}  {"Output Stage":^15}')
    print(f'{"===============":15}  {"===============":15}  {"===============":15}  {"===============":15}')
    for k,v in res.items():
        print(f'{nhd_feature_id:<15}  {k:<15}  {CMS:<15}  {v:<15}')
