#!/usr/bin/env python3

"""
Purpose: Utility functions for generating FIM maps using HAND published
         by NOAA OWP.
Authors: Tony Castronova <acastronova@cuahsi.org>
Last Modified: March 13, 2024
"""


import numpy
import xarray
import rioxarray
import geopandas
import fim_utils as fim
from pathlib import Path
import matplotlib.pyplot as plt
from geocube.api.core import make_geocube

def get_hand_object(path):
    return rioxarray.open_rasterio(path, masked=True).squeeze().drop('band').to_dataset(name='hand')

def generate_fim_grid(nhd_feature_id: int,
                      cms: float) -> xarray.Dataset:
    
    # Collect stage for all hydroids in this NHD+ feature. 
    # Use the input cms to interpolate river stage from a rating curve.
    stage_dict = fim.get_stage_for_all_hydroids_in_reach(nhd_feature_id, cms)
    
    # Load the precomputed HAND raster.
    xds = get_hand_object(Path('./rem_zeroed_masked_0.tif'))
    
    # Make a copy of the 'hand' variable to the 'stage' variable.
    # This will be used to compute the flood inundation map.
    xds['stage'] = xds.hand.copy(deep=True)

    # Read watershed geometries and set stage values from the `stage_dict` defined above.
    # Remove all other geometries.
    geodf = geopandas.read_file(Path('./gw_catchments_reaches_filtered_addedAttributes_crosswalked_0.gpkg'))

    for hydroid, stage in stage_dict.items():
        geodf.loc[geodf.HydroID==hydroid, 'stage'] = stage

    geodf_filtered = geodf[geodf.stage.notnull()]

    # Add these geometries to the HAND and Stage DataSet using a GeoCube.
    # Create variables for stage and hydroID to the HAND raster, drop 
    # everything except the HydroIDs that we're interested in.
    
    out_grid = make_geocube(
        vector_data=geodf,
        measurements=['HydroID'],
        like=xds # ensure the data are on the same grid
    ) 
    xds = xds.assign_coords( hydroid = (['y', 'x'], out_grid.HydroID.data) )
    xds = xds.where(xds.hydroid.isin(geodf_filtered.HydroID), drop=True)

    # Update the stage values in the DataSet where specific hydroid's exist. 
    for idx, row in geodf_filtered.iterrows():
        xds['stage'] = xarray.where(xds.hydroid == row.HydroID, row.stage, xds.stage)

    # Compute FIM by subtracting `hand` from `stage`.
    # Everything that is negative should be set to zero.
    xds['fim'] = xds.stage - xds.hand
    xds['fim'] = xarray.where(xds.fim >= 0.00001, xds.fim, numpy.nan)
    
    return xds

if __name__ == '__main__':
    
    nhd_feature_id=946010122
    cms = 400
    generate_fim_grid(nhd_feature_id, cms)