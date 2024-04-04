

import pytest
from nwm_utils import fim


class TestXFim:
    nhd_feature_id = 946010122
    st_date = '2000-03-01'
    end_date = '2020-06-30'

    def test_generate_fim_grid(self):
        cms = 100.
        ds = fim.xfim.generate_fim_grid(self.nhd_feature_id,
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
        ds = fim.xfim.generate_fim_through_time(self.nhd_feature_id,
                                               times,
                                               flows)

        assert (ds.sizes['x'], ds.sizes['y'], ds.sizes['time']) == (345, 612, 5)
        assert ds.fim_extent.sum().item() == 148362.21875
        


