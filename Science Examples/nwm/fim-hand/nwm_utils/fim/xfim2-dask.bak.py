import dask
import pandas
import numpy as np
import xarray as xr
import geopandas as gpd
from pathlib import Path
from typing import List, Dict
from datetime import datetime
from dask.delayed import delayed
from dask.distributed import Client, print
from geocube.api.core import make_geocube
from collections import OrderedDict
import utils
from utils import get_hand_object, get_stage_for_all_hydroids_in_reach

@dask.delayed
def get_stage_from_rating_curve(nhd_feature_id: int, flow: float) -> Dict[int, float]:
    """
    Collects stage for all hydroids within the given NHD feature.
    Stage is computed by interpolating known rating curves using
    the provided `flow`.

    Parameters:
    nhd_feature_id: int
        NHD+ feature identifier for which to compute FIM.
    flow: float
        Streamflow (in cubic meters per second) used to derive river stage.

    Returns:
    Dict[int, float] | delayed
        Dictionary containing one or more NWM hydro identifiers and
        their corresponding stages.
    """
    return delayed(get_stage_for_all_hydroids_in_reach)(nhd_feature_id, flow)

@dask.delayed
def set_stage_value(ds: xr.Dataset, stage_dict: Dict[int, float] | delayed) -> xr.DataArray:
    """
    Sets stage in Xarray Dataset with or without Dask.

    Parameters:
    ds: xr.Dataset
        The input dataset.
    stage_dict: Dict[int, float] | delayed
        Dictionary containing NWM hydro identifiers and their corresponding stages.

    Returns:
    xr.DataArray
        Array containing stage values.
    """
    a = np.empty(ds.hand.shape)
    a[:] = np.nan

    hydroids = stage_dict.keys()
    for hydroid in hydroids:
        a[np.where(ds.hydroid == hydroid)] = stage_dict[hydroid]

    da = xr.DataArray(data=a, dims=['y', 'x'], coords=dict(x=ds.x, y=ds.y))
    return da

@dask.delayed
def build_stage_data_array(ds: xr.Dataset, reaches: List[utils.Reach]) -> xr.DataArray:
    """
    This sets the stages for all times at each reach.

    Parameters:
    ds: xr.Dataset
        The input dataset.
    reaches: List[utils.Reach]
        List of Reach objects.

    Returns:
    xr.DataArray
        Array containing stage data for all times at each reach.
    """
    dat = OrderedDict()
    for reach in reaches:
        rcs = []
        hydro_ids = []
        for hydro_id, rc in reach.get_rating_curve().items():
            hydro_ids.append(hydro_id)
            rcs.append(rc)

        for time, flow in zip(reach.times(), reach.flows()):
            if time not in dat:
                dat[time] = {}
            for i in range(len(hydro_ids)):
                stage = rcs[i].get_stage_for_flow(flow)
                dat[time][hydro_ids[i]] = stage
    
    futures = [set_stage_value(ds, dat[time]) for time in dat.keys()]
    das = dask.compute(*futures)
    da = xr.concat(das, pandas.DatetimeIndex(dat.keys(), name='time'))
    
    return da

@dask.delayed
def set_stage_for_hydroids(ds: xr.Dataset, reaches: List[utils.Reach]) -> xr.Dataset:
    """
    Set stage values in the dataset for all reaches.

    Parameters:
    ds: xr.Dataset
        The input dataset.
    reaches: List[utils.Reach]
        List of Reach objects.

    Returns:
    xr.Dataset
        The modified dataset with stage values.
    """
    for reach in reaches:
        dat = OrderedDict()
        for time, flow in zip(reach.times(), reach.flows()):
            if time not in dat:
                dat[time] = {}
            for hydro_id, rc in reach.get_rating_curve().items():
                stage = rc.get_stage_for_flow(flow)
                dat[time][hydro_id] = stage
        print('HELLO')
        da = build_stage_data_array(ds, [reach])
        #ds['stage'] = da
    return ds

@dask.delayed
def generate_fim(ds: xr.Dataset, times: List[datetime]) -> xr.Dataset:
    """
    Generates FIM grids for a range of times.

    Parameters:
    ds: xr.Dataset
        The input dataset.
    times: List[datetime]
        List of times for which FIM will be created.

    Returns:
    xr.Dataset
        Dataset containing FIM grids.
    """
    fims = []
    masks = []
    for time in times:
        da_fim, da_fim_mask = compute_fim(ds.sel(time=time))
        fims.append(da_fim)
        masks.append(da_fim_mask)

    fim_array = xr.concat(fims, dim='time')
    mask_array = xr.concat(masks, dim='time')

    ds['fim'] = fim_array
    ds['fim_mask'] = mask_array
    return ds

@dask.delayed
def compute_fim(ds: xr.Dataset) -> xr.DataArray:
    """
    Computes FIM for a single time step.

    Parameters:
    ds: xr.Dataset
        The input dataset for a single time step.

    Returns:
    xr.DataArray
        DataArray containing FIM for the given time step.
    """
    da_fim = ds['stage'] - ds['hand']
    da_fim = xr.where(da_fim >= 0.00001, da_fim, np.nan)
    da_fim_mask = xr.where(da_fim >= 0.00001, 1, da_fim)
    da_fim_mask = xr.where(da_fim < 0.00001, 0, da_fim_mask)
    return da_fim, da_fim_mask

if __name__ == '__main__':
    client = Client()

    reach_ids = [946010122, 946010121, 946010120]
    st_date = '2020-03-01'
    end_date = '2020-03-10'

    xds = get_hand_object(Path("/Users/castro/Documents/work/notebooks/nwm/fim-hand/rem_zeroed_masked_0.tif"))

    gdf = gpd.read_file(Path("/Users/castro/Documents/work/notebooks/nwm/fim-hand/gw_catchments_reaches_filtered_addedAttributes_crosswalked_0.gpkg"))
    subset_gdf = gdf.loc[gdf.feature_id.isin(reach_ids)]
    xds = xds.rio.clip(subset_gdf.geometry.values, gdf.crs, drop=True, invert=False, from_disk=True, all_touched=True)
    gdf = gdf.loc[gdf.feature_id.isin(reach_ids)]

    out_grid = make_geocube(vector_data=gdf, measurements=["HydroID"], like=xds)
    xds = xds.assign_coords(hydroid=(["y", "x"], out_grid.HydroID.data))
    xds = xds.where(xds.hydroid.isin(gdf.HydroID), drop=True)

    dates = pandas.date_range(st_date, end_date, freq='h')
    flows = [np.random.uniform(40, 300) for _ in range(len(dates))]

    hydro_df = pandas.read_csv('/Users/castro/Documents/work/notebooks/nwm/fim-hand/hydroTable_0.csv', usecols=["HydroID", "NextDownID", "feature_id", "stage", "discharge_cms"])

    reaches = [utils.Reach(reach_id, flows, dates) for reach_id in reach_ids]
    for reach in reaches:
        reach.set_rating_curve_from_hydrotable(hydro_df)

    xds2 = set_stage_for_hydroids(xds, reaches).compute()
    #xds2 = generate_fim(xds2, dates).compute()

    print(xds2)
