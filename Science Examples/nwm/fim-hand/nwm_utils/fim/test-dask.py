

import pytest

import dask
import numpy
import xarray
import geopandas
from pathlib import Path
from nwm_utils import fim
from dask.delayed import Delayed
from geocube.api.core import make_geocube


class FIMTestData():
    def __init__(self,
                 nhd_feature_id,
                 st_date,
                 end_date,
                 hand_raster,
                 xds):
        self.nhd_feature_id = nhd_feature_id
        self.st_date = st_date
        self.end_date = end_date
        self.hand_raster = hand_raster
        self.xds = xds


@pytest.fixture(scope='session', autouse=True)
def initialize():

    nhd_feature_id = 946010122
    fpath = Path(__file__).resolve().parent
    hand_raster = fpath/Path('../tests/test-data/rem_zeroed_masked_0.tif')

    xds = fim.xfim.get_hand_object(hand_raster)

    gdf = geopandas.read_file(
        fpath/Path("../tests/test-data/gw_catchments_reaches_filtered_addedAttributes_crosswalked_0.gpkg"))

    subset_gdf = gdf.loc[gdf.feature_id == nhd_feature_id]
    xds = xds.rio.clip(subset_gdf.geometry.values,
                       gdf.crs,
                       drop=True,
                       invert=False,
                       from_disk=True,
                       all_touched=True)

    gdf = gdf[gdf.feature_id == nhd_feature_id]

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

    pytest.nhd_feature_id = nhd_feature_id
    pytest.st_date = '2000-03-01'
    pytest.end_date = '2020-06-30'
    pytest.hand_raster = hand_raster
    pytest.xds = xds


class TestXFim:

    def test_set_stage_in_dataset(self):
        cms = 100
        from datetime import datetime
        times = [datetime(2000,1,1),
                 datetime(2000,1,2),
                 datetime(2000,1,3),
                 datetime(2000,1,4),
                 datetime(2000,1,5)]

        flows = [cms] * len(times)
        reach_ids = [pytest.nhd_feature_id]

        # generate some random flow values for each of these hydroids
        flow_obj = {}
        for i in range(0, len(times)):
            flow_obj[times[i]] = [flows[i]]

        ds = fim.xfim.set_stage_for_hydroids(pytest.xds,
                                             reach_ids,
                                             flow_obj)
        del client
        
        import pandas
        compare_dates = [pandas.to_datetime(t) for t in times]
        assert list(ds.time.values) == compare_dates

    def test_generate_fim(self):
        cms = 100
        from datetime import datetime
        times = [datetime(2000,1,1),
                 datetime(2000,1,2),
                 datetime(2000,1,3),
                 datetime(2000,1,4),
                 datetime(2000,1,5)]

        flows = [cms] * len(times)
        reach_ids = [pytest.nhd_feature_id]

        # generate some random flow values for each of these hydroids
        flow_obj = {}
        for i in range(0, len(times)):
            flow_obj[times[i]] = [flows[i]]

       
        # set stage for hydroids
        ds = fim.xfim.set_stage_for_hydroids(pytest.xds, reach_ids, flow_obj)

        # scatter the data ahead of time so we can pass a pointer to the 
        # scheduler instead of the entire data object
        scattered_xds = client.scatter(ds, broadcast=True)
 
        ds = fim.xfim.generate_fim(scattered_xds, times, parallel=True)

        del client
        import pandas

        compare_dates = [pandas.to_datetime(t) for t in times]
        assert list(ds.time.values) == compare_dates
        assert int(ds.fim.sum().item()) == 148362



