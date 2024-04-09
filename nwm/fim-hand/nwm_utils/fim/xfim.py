#!/usr/bin/env python3

"""
Purpose: Utility functions for generating FIM maps using HAND published
         by NOAA OWP.
Authors: Tony Castronova <acastronova@cuahsi.org>
Last Modified: April 09, 2024
"""

import dask
import numpy
import xarray
import pandas
import geopandas
import fim_utils as fim
from pathlib import Path
from typing import List, Dict
from datetime import datetime
from dask.delayed import Delayed
import dask.distributed as distributed
from geocube.api.core import make_geocube

from .utils import get_hand_object, get_stage_for_all_hydroids_in_reach


__all__ = ['generate_fim_grid']


def get_stage_from_rating_curve(nhd_feature_id: int,
                                flow: float,
                                delayed: bool = False) -> Dict [numpy.int64, float] | Delayed:
    """
    Collects stage for all hydroids within the given nhd feature.
    Stage is computed by iterpolating known rating curves using
    the provide `flow`.

    Parameters
    ==========
    nhd_feature_id: int
        NHD+ feature identifier for which to compute FIM.
    cms: float
        Streamflow (in cubic meters per second) used to derive river stage.
    delayed: bool = False
        Flag to indicate if this function should use Dask for parallelization.

    Returns
    =======
    Dict [int, float] | dask.delayed.Delayed
        Dictionary containing one or more NWM hydro identifiers and
        their corresponding stages.

    """
    if delayed:
        return __dask_get_stage_from_rating_curve(nhd_feature_id, flow)
    return get_stage_for_all_hydroids_in_reach(nhd_feature_id, flow)


@dask.delayed
def __dask_get_stage_from_rating_curve(nhd_feature_id: int,
                                flow: float) -> Delayed:
    """
    Collects stage for all hydroids within the given nhd feature.
    Stage is computed by iterpolating known rating curves using
    the provide `flow`.
    """
    
    return get_stage_for_all_hydroids_in_reach(nhd_feature_id, flow)

def set_stage_value(ds: xarray.Dataset,
                    stage_dict: Dict [numpy.int64, float] | Delayed,
                    delayed=False):
    """
    Sets stage in Xarray Dataset with out without dask.
    """
    
    if delayed:
        return __dask_set_stage_value(ds, stage_dict)
    return __set_stage_value(ds, stage_dict)

def __set_stage_value(ds: xarray.Dataset,
                    stage_dict: Dict [numpy.int64, float] | Delayed,
                    delayed=False):
    """
    Sets stage in Xarray Dataset.
    """
    
    if delayed:
        return __dask_set_stage_value(ds, stage_dict)

    # create an array to store stage values for the current timestep
    a = numpy.empty(ds.hand.shape)
    a[:] = numpy.nan
    
    # Update the stage values in the array where specific hydroid's exist.
    hydroids = stage_dict.keys()
    for hydroid in hydroids:
        a[numpy.where(ds.hydroid == hydroid)] = stage_dict[hydroid]
    
    da = xarray.DataArray(data=a, dims=['y','x'],
                          coords = dict(
                             x=ds.x,
                             y=ds.y
                         ))
    return da

@dask.delayed
def __dask_set_stage_value(ds: xarray.Dataset,
                           stage_dict: Dict [numpy.int64, float] | Delayed):
    """
    Returns a future for setting stage in Xarray Dataset.
    """
   
    return __set_stage_value(ds, stage_dict)


def build_stage_data_array(ds, reach_ids, flows,
                           delayed=True):
    """
    This sets the stages for all times at each reach.
    """

    if delayed:
        return __dask_build_stage_data_array(ds, reach_ids, flows)

    return __build_stage_data_array(ds, reach_ids, flows)

def __build_stage_data_array(ds, reach_ids, flows):
    """
    This should set the stages for all times at each reach.
    """
    return Exception('Non-delayed functionality not Implemented')

@dask.delayed
def __dask_build_stage_data_array(ds, reach_ids, flows):
    """
    This should set the stages for all times at each reach.
    """
    
    
    futures = []
    for time in flows.keys():
        # the array of flows for each reach at the current timestep
        reach_flows = flows[time]
        
        # compute stage for each reach in the current timestep
        stage_dict = {}
        
        for i in range(0, len(reach_ids)): 
            nhd_feature_id = reach_ids[i]
            flow = reach_flows[i]
            futures.append(get_stage_from_rating_curve(nhd_feature_id, flow))
            
    stages = dask.compute(futures)[0]

    i = 0
    futures = []
    for time in flows.keys():
        stage_dict = stages[i]
        futures.append(set_stage_value(ds, stage_dict))
        i += 1
    das = dask.compute(futures)[0]

    da = xarray.concat(das,
                       pandas.DatetimeIndex(flows.keys(), name='time'))

    return da

def set_stage_for_hydroids(ds, reach_ids, flows, parallel=True):

    if parallel:
        client = distributed.client._get_global_client() or distributed.Client()

        # scatter the data ahead of time so we can pass a pointer to the 
        # scheduler instead of the entire data object
        scattered_ds = client.scatter(ds, broadcast=True)
        
        res = dask.compute(build_stage_data_array(scattered_ds,
                                                  reach_ids,
                                                  flows,delayed=True,
                                                  ))[0]

        # set result as variable in the input dataset
        ds['stage'] = res
        return ds
    else:
        raise Exception('Non-parallel functionality not implemented yet')

def generate_fim(ds, times, parallel=True):
    
    if parallel:
        
        return __dask_generate_fim(ds, times, parallel=True)
    else:
        return Exception("Non-parallel is not implemented")



def __dask_generate_fim(xds, times, parallel=True):


#    client = distributed.client._get_global_client() or distributed.Client()
    
#    # scatter the data ahead of time so we can pass a pointer to the 
#    # scheduler instead of the entire data object
#    scattered_xds = client.scatter(xds, broadcast=True)

    futures = []
    for time in times:
        future = compute_fim(xds,
                             time)
        futures.append(future)
    res = dask.compute(futures)

    fims, masks = zip(*res[0])
    
    fims = xarray.concat(fims, dim='time')
    masks = xarray.concat(masks, dim='time')

    # convert lazy collection into dask collection this
    # is necessary if xds has been scattered.
    if isinstance(xds, dask.distributed.client.Future):
        xds = xds.result()

    xds['fim'] = fims
    xds['fim_mask'] = masks
    
    return xds



@dask.delayed
def compute_fim(xds, time):
    
    da_fim = (xds.sel(time = time).stage - xds.hand)
    da_fim = xarray.where(da_fim >= 0.00001, da_fim, numpy.nan)
    da_fim_mask = xarray.where(da_fim >= 0.00001, 1, da_fim)
    da_fim_mask = xarray.where(da_fim < 0.00001, 0, da_fim_mask)
    
    return da_fim, da_fim_mask

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
    xds = get_hand_object(Path("./rem_zeroed_masked_0.tif"))

    # Make a copy of the 'hand' variable to the 'stage' variable.
    # This will be used to compute the flood inundation map.
    xds["stage"] = xds.hand.copy(deep=True)

    # Read watershed geometries and set stage values from the `stage_dict` defined above.
    # Remove all other geometries.
    geodf = geopandas.read_file(
        Path("./gw_catchments_reaches_filtered_addedAttributes_crosswalked_0.gpkg")
    )

    for hydroid, stage in stage_dict.items():
        geodf.loc[geodf.HydroID == hydroid, "stage"] = stage

    geodf_filtered = geodf[geodf.stage.notnull()]

    # Add these geometries to the HAND and Stage DataSet using a GeoCube.
    # Create variables for stage and hydroID to the HAND raster, drop
    # everything except the HydroIDs that we're interested in.

    out_grid = make_geocube(
        vector_data=geodf,
        measurements=["HydroID"],
        like=xds,  # ensure the data are on the same grid
    )
    xds = xds.assign_coords(hydroid=(["y", "x"], out_grid.HydroID.data))
    xds = xds.where(xds.hydroid.isin(geodf_filtered.HydroID), drop=True)

    # Update the stage values in the DataSet where specific hydroid's exist.
    for _, row in geodf_filtered.iterrows():
        xds["stage"] = xarray.where(xds.hydroid == row.HydroID, row.stage, xds.stage)

    # Compute FIM by subtracting `hand` from `stage`.
    # Everything that is negative should be set to zero.
    xds["fim"] = xds.stage - xds.hand
    xds["fim_extent"] = xarray.where(xds.fim >= 0.00001, xds.fim, numpy.nan)
    xds["fim_mask"] = xarray.where(xds.fim >= 0.00001, 1, xds.fim)
    xds["fim_mask"] = xarray.where(xds.fim < 0.00001, 0, xds.fim_mask)

    return xds
    

def generate_fim_through_time(nhd_feature_id: int,
                              times: List[datetime],
                              flows: List[float]) -> xarray.Dataset:
    """
    Computes FIM grids for a range of times at a single NHD+ reach.

    Parameters
    ==========
    nhd_feature_id: int
        NHD+ feature identifier for which to compute FIM.
    times: List[datetime.datetime]
        List of times for which FIM will be created
    flows: List[float]
        List of streamflows (in cubic meters per second) used to derive river stage.

    Returns
    =======
    xarray.Dataset
        Dataset containing HAND, FIM, and Stage variables through time
    """
    
    hand = None
    datasets = []
    for i, dat in enumerate(zip(times, flows)):
        time = dat[0]
        cms = dat[1]
        print(f'Computing FIM: NHD+ {nhd_feature_id} at {time}')
        
        # compute fim
        ds = generate_fim_grid(nhd_feature_id, cms)
        
        # save hand grid for the first iteration only.
        # this will be re-added after concatenating
        # all of the results
        if i == 0:
            hand = ds.hand.copy(deep=True)
        
        # drop the hand variable because we don't want
        # this to be duplicated through time
        ds.drop_vars(['hand'])
        
        # save results for concatenation later
        datasets.append(ds)
    
    concat_ds = xarray.concat(datasets,
                              pandas.DatetimeIndex(times, name='time'))
    concat_ds['hand'] = hand
    
    return concat_ds
