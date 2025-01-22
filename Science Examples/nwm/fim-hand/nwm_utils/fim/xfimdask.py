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
from pathlib import Path
from typing import List, Dict
from datetime import datetime
from dask.delayed import Delayed
import dask.distributed as distributed
from geocube.api.core import make_geocube

from collections import OrderedDict

import utils
from utils import get_hand_object, get_stage_for_all_hydroids_in_reach


__all__ = ['generate_fim_grid']


def get_stage_from_rating_curve(nhd_feature_id: int,
                                flow: float) -> Dict [numpy.int64, float] | Delayed:
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

    return get_stage_for_all_hydroids_in_reach(nhd_feature_id, flow)

def set_stage_value(ds: xarray.Dataset,
                    stage_dict: Dict [numpy.int64, float] | Delayed):
    """
    Sets stage in Xarray Dataset with out without dask.
    """

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


def build_stage_data_array(ds, reaches):
    """
    This sets the stages for all times at each reach.
    """

    # loop through time keys
    das = []
    for time in reaches.keys():
        das.append(set_stage_value(ds, reaches[time]))
    da = xarray.concat(das,
                       pandas.DatetimeIndex(reaches.keys(), name='time'))
    return da

def set_stage_for_hydroids(ds, reaches):


    # collect hydro reaches and their rating curves
    print('preparing data')
    dat = OrderedDict()
    for reach in reaches:
        rcs = []
        hydro_ids = []
        for hydro_id, rc in reach.get_rating_curve().items():
            hydro_ids.append(hydro_id)
            rcs.append(rc)

        # loop through reach time(s) and flow(s), compute stage
        for _, data in enumerate(zip(reach.times(), reach.flows())):
            time = data[0]
            flow = data [1]
        
            if time not in dat:
                dat[time] = {}
                
            for i in range(0, len(hydro_ids)):
                stage = rcs[i].get_stage_for_flow(flow)
                dat[time][hydro_ids[i]] = stage

        print('building stage array') 
        da = build_stage_data_array(ds, dat)
        # set result as variable in the input dataset
        ds['stage'] = da

    return ds


def generate_fim(ds, times):
    
    res = []
    for time in times:
        res.append(compute_fim(ds,
                             time))

    fims, masks = zip(*res)
    
    fims = xarray.concat(fims, dim='time')
    masks = xarray.concat(masks, dim='time')

    xds['fim'] = fims
    xds['fim_mask'] = masks
    
    return xds

def compute_fim(ds, time):
    
    da_fim = (ds.sel(time = time).stage - ds.hand)
    da_fim = xarray.where(da_fim >= 0.00001, da_fim, numpy.nan)
    da_fim_mask = xarray.where(da_fim >= 0.00001, 1, da_fim)
    da_fim_mask = xarray.where(da_fim < 0.00001, 0, da_fim_mask)
    
    return da_fim, da_fim_mask

if __name__ == '__main__':


    # use a try accept loop so we only instantiate the client
    # if it doesn't already exist.
    from dask.distributed import Client
    try:
        print(client.dashboard_link)
    except:    
        # The client should be customized to your workstation resources.
        client = Client(n_workers=2, memory_limit='4GB') # per worker
        print(client.dashboard_link)
    
    import time
    st = time.time()
    
    reach_ids = [946010122, 946010121, 946010120]
    st_date = '2020-03-01'
    end_date = '2020-03-10'
    
    # reduce the size of the data we're working with
    print('preparing data')
    xds = get_hand_object(Path("/Users/castro/Documents/work/notebooks/nwm/fim-hand/rem_zeroed_masked_0.tif"))
    
    gdf = geopandas.read_file(
            Path("/Users/castro/Documents/work/notebooks/nwm/fim-hand/gw_catchments_reaches_filtered_addedAttributes_crosswalked_0.gpkg")
        )
    subset_gdf = gdf.loc[gdf.feature_id.isin(reach_ids)]
    xds = xds.rio.clip(subset_gdf.geometry.values,
                     gdf.crs,
                     drop=True,
                     invert=False,
                     from_disk=True,
                     all_touched=True
                    )
    
    gdf = gdf.loc[gdf.feature_id.isin(reach_ids)]
    
    # Add these geometries to the HAND and Stage DataSet using a GeoCube.
    # Create variables for stage and hydroID to the HAND raster, drop
    # everything except the HydroIDs that we're interested in.
    out_grid = make_geocube(
        vector_data=gdf,
        measurements=["HydroID"],
        like=xds,  # ensure the data are on the same grid
    )
    xds = xds.assign_coords(hydroid=(["y", "x"], out_grid.HydroID.data))
    xds = xds.where(xds.hydroid.isin(gdf.HydroID), drop=True)


    # make up some data
    import random
    dates = pandas.date_range(st_date, end_date, freq='h')
    flows = [random.uniform(40, 300) for _ in range(len(dates))]

    hydro_df = pandas.read_csv(
    '/Users/castro/Documents/work/notebooks/nwm/fim-hand/hydroTable_0.csv',
    usecols=["HydroID", "NextDownID", "feature_id", "stage", "discharge_cms"],
    )

    print('setting stage in class')
    reaches = []
    for reach_id in reach_ids:
        print(f'{reach_id}')
        r = utils.Reach(reach_id, flows, dates)
        r.set_rating_curve_from_hydrotable(hydro_df)
        reaches.append(r)
        print('---')

    print('setting stage in xarray')
    xds2 = set_stage_for_hydroids(xds, reaches)

    print('generate fim')
    xds2 = generate_fim(xds2, dates)

    print(xds2)
    print(f'Elapsed time = {time.time() - st}')

