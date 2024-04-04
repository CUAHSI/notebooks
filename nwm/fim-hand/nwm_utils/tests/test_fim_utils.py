
import pdb

import pytest
from pathlib import Path
from nwm_utils import fim

class TestFimUtils:

    fpath = Path(__file__).resolve().parent
    hand_raster = fpath/Path('test-data/rem_zeroed_masked_0.tif')
    nhd_feature_id = 946010122
    geo_lon = -111.0538022
    geo_lat = 42.21104338


    def test_get_hand_object(self):
        ds = fim.utils.get_hand_object(self.hand_raster)
        assert (ds.sizes['x'], ds.sizes['y']) == (6285, 8696)
        assert ds.hand.sum().item() == 3453847040.0

#    def test_get_nwm_id(self):
#        import xarray
#        import fsspec
#
#        ds = conus_bucket_url = 's3://noaa-nwm-retrospective-3-0-pds/CONUS/zarr/chrtout.zarr'
#        ds = xarray.open_zarr(fsspec.get_mapper(
#                                conus_bucket_url,
#                                anon=True))
#
#        nwm_id = fim.utils.get_nwm_id(self.geo_lat, self.geo_lon, ds)
#        pdb.set_trace()


    def test_get_single_stage(self):
        flow = 100
        dat = fim.utils.get_stage_for_nhd_reach(self.nhd_feature_id, flow)

        assert 27670095 in dat.keys()
        assert dat[27670095] == 3.903178622396883

    def test_get_all_stages_for_reach(self):
        flow = 100
        dat = fim.utils.get_stage_for_all_hydroids_in_reach(self.nhd_feature_id, flow)

        assert len(dat.keys()) == 7
        assert sum(dat.values()) == 32.36753983826069

