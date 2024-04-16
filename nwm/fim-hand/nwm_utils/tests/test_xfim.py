

import pytest

import dask
import numpy
import xarray
import pandas
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
                 xds,
                 hydrodf):
        self.nhd_feature_id = nhd_feature_id
        self.st_date = st_date
        self.end_date = end_date
        self.hand_raster = hand_raster
        self.xds = xds
        self.hydrodf = hydrodf


@pytest.fixture(scope='session', autouse=True)
def initialize():

    nhd_feature_id = 946010122
    fpath = Path(__file__).resolve().parent
    hand_raster = fpath/Path('test-data/rem_zeroed_masked_0.tif')

    xds = fim.xfim.get_hand_object(hand_raster)

    gdf = geopandas.read_file(
        fpath/Path("test-data/gw_catchments_reaches_filtered_addedAttributes_crosswalked_0.gpkg"))
    
    hydrodf = pandas.read_csv(fpath/Path('test-data/hydroTable_0.csv'),
                              usecols=["HydroID",
                                       "NextDownID",
                                       "feature_id",
                                       "stage",
                                       "discharge_cms"])


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
    pytest.hydrodf = hydrodf


class TestXFim:


    def test_generate_fim_grid(self):
        cms = 100.
        ds = fim.xfim.generate_fim_grid(pytest.nhd_feature_id,
                                        cms)
       
        assert (ds.sizes['x'], ds.sizes['y']) == (345, 612)
        assert ds.fim_extent.sum().item() == 29672.453125


    def test_generate_fim_through_time(self):
        cms = 100.

        from datetime import datetime
        times = [datetime(2000,1,1),
                 datetime(2000,1,2),
                 datetime(2000,1,3),
                 datetime(2000,1,4),
                 datetime(2000,1,5)]

        flows = [cms] * len(times)
        ds = fim.xfim.generate_fim_through_time(pytest.nhd_feature_id,
                                               times,
                                               flows)

        assert (ds.sizes['x'], ds.sizes['y'], ds.sizes['time']) == (345, 612, 5)
        assert ds.fim_extent.sum().item() == 148362.21875

    def test_get_stage(self):
        cms = 100.
        stage_dict = fim.xfim.get_stage_from_rating_curve(pytest.nhd_feature_id,
                                                          flow=cms)
        assert list(set(map(type, stage_dict.keys())))[0] == numpy.int64
        assert list(set(map(type, stage_dict.values())))[0] == float
        assert len(stage_dict.values()) == 7


        stage_future = fim.xfim.get_stage_from_rating_curve(pytest.nhd_feature_id,
                                                            flow=cms,
                                                            delayed=True)
        assert type(stage_future) == Delayed 

        stage_dict = dask.compute(stage_future)[0]
        assert list(set(map(type, stage_dict.keys())))[0] == numpy.int64
        assert list(set(map(type, stage_dict.values())))[0] == float
        assert len(stage_dict.values()) == 7


    def test_set_stage(self):
        cms = 100.
        stage_dict = fim.xfim.get_stage_from_rating_curve(pytest.nhd_feature_id,
                                                          flow=cms)
        da = fim.xfim.set_stage_value(pytest.xds, stage_dict)

        assert type(da) == xarray.DataArray

        delayed = fim.xfim.set_stage_value(pytest.xds, stage_dict, delayed=True)
        assert type(delayed) == Delayed

        da2 = dask.compute(delayed)[0]
        assert da.broadcast_equals(da2)


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

        reaches = []        
        for reach_id in reach_ids:
            r = fim.utils.Reach(reach_id, flows, times)
            r.set_rating_curve_from_hydrotable(pytest.hydrodf)
            reaches.append(r)

#        from dask.distributed import Client
#        client = Client(n_workers=2, memory_limit='4GB') # per worker

        # scatter the data ahead of time so we can pass a pointer to the 
        # scheduler instead of the entire data object
#        scattered_ds = client.scatter(pytest.xds, broadcast=True)
        ds = fim.xfim.set_stage_for_hydroids(pytest.xds,
                                             reaches,
                                             parallel=False)
 #       del client
        
        compare_dates = [pandas.to_datetime(t) for t in times]
        assert list(ds.time.values) == compare_dates

        # TODO: test that stage values are set correctly
        

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

        reaches = []        
        for reach_id in reach_ids:
            r = fim.utils.Reach(reach_id, flows, times)
            r.set_rating_curve_from_hydrotable(pytest.hydrodf)
            reaches.append(r)

#        from dask.distributed import Client
#        client = Client(n_workers=2, memory_limit='4GB') # per worker
       
        # set stage for hydroids
        ds = fim.xfim.set_stage_for_hydroids(pytest.xds, reaches, parallel=False)
        import pdb; pdb.set_trace()
        # scatter the data ahead of time so we can pass a pointer to the 
        # scheduler instead of the entire data object
#        scattered_xds = client.scatter(ds, broadcast=True)
 
        ds = fim.xfim.generate_fim(ds, times, parallel=False)

 #       del client
        compare_dates = [pandas.to_datetime(t) for t in times]
        assert list(ds.time.values) == compare_dates
        assert int(ds.fim.sum().item()) == 148362

    def test_compute_fim(self):
        pass



