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
from geocube.api.core import make_geocube

def get_hand_object(path: Path) -> xarray.Dataset:
    """
    Converts NOAA OWP HAND tif into a rioxarray Dataset.

    Parameters
    ==========
    path: pathlib.Path
        Path-like object to the NOAA OWP HAND tif.

    Returns
    =======
    xarray.Dataset
        Dataset containing HAND as a variable.

    """
    return rioxarray.open_rasterio(path, masked=True) \
            .squeeze() \
            .drop_vars('band'). \
            to_dataset(name='hand')

def generate_fim_grid(nhd_feature_id: int,
                      cms: float) -> xarray.Dataset:
    """
    Computes FIM grid for a given NHD+ feature ID and flow rate.

    Parameters
    ==========
    nhd_feature_id: int
        NHD+ feature identifier for which to compute FIM.
    cms: float
        Streamflow (in cubic meters per second) used to derive river stage.

    Returns
    =======
    xarray.Dataset
        Dataset containing HAND, FIM, and Stage variables.
    """
    
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
    for _, row in geodf_filtered.iterrows():
        xds['stage'] = xarray.where(xds.hydroid == row.HydroID, row.stage, xds.stage)

    # Compute FIM by subtracting `hand` from `stage`.
    # Everything that is negative should be set to zero.
    xds['fim'] = xds.stage - xds.hand
    xds['fim_extent'] = xarray.where(xds.fim >= 0.00001, xds.fim, numpy.nan)
    xds['fim_mask'] = xarray.where(xds.fim >= 0.00001, 1, xds.fim)
    xds['fim_mask'] = xarray.where(xds.fim < 0.00001, 0, xds.fim_mask)

    return xds

if __name__ == '__main__':
    
    nhd_feature_id=946010122
    cms = 400
    generate_fim_grid(nhd_feature_id, cms)
