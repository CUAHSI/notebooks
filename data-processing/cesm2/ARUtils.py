# /usr/bin/env python3


import numpy as np
import xarray as xr


# Global params
AR_THRESHOLD = 250  # minimum threshold for AR detection (250 kg/m/s)

# Initialize AR Category lookup table
# AR[IVT/10, Duration (hours)]
AR = np.zeros((200, 144)) + 5
# cat 4
AR[125:150, :24] = 4
AR[100:125, 24:48] = 4
AR[75:100, 48:72] = 4
AR[50:75, 72:96] = 4
AR[25:50, 96:120] = 4
# cat 3
AR[100:125, :24] = 3
AR[75:100, 24:48] = 3
AR[50:75, 48:72] = 3
AR[25:50, 72:96] = 3
# cat 2
AR[75:100, :24] = 2
AR[50:75, 24:48] = 2
AR[25:50, 48:72] = 2
# cat 1
AR[50:75, :24] = 1
AR[25:50, 24:48] = 1
# cat 0
AR[:25, :] = 0
AR[25:50, :24] = 0


def create_ar_index_for_grid_cell(ds):
    """
    This function computes an AR index and AR duration for a given xarray array.

    Parameters:
    ===========
    ds: xarray Dataset containg IVT

    Returns:
    ========
    xarray_array with AR index and AR duration variables
    """

    # isolate the days in which the minimum AR conditions are met.
    ds_thr = xr.where((ds.IVT >= AR_THRESHOLD), 1, 0)

    # Flatten the array
    flattened = ds_thr.values
    flattened = np.insert(flattened, 0, 0, axis=0)
    flattened = np.append(flattened, 0)

    ar_events = np.zeros(flattened.shape)
    ar_duration = np.zeros(flattened.shape)

    # Initialize variables
    event = 0
    idx_start = 0

    # Initialize variables
    current_sequence = 0

    for idx in range(1, len(flattened)):

        # we hit a value of 1
        if flattened[idx] == 1:
            current_sequence += 1

            if flattened[idx - 1] == 0:
                # this is the first value in the set
                idx_start = idx
                event += 1

        # we hit a value of 0
        else:
            if current_sequence > 0:
                ar_duration[idx_start:idx] = current_sequence
            current_sequence = 0
            ar_events[idx_start:idx] = event

    # return all but the first zero that we added
    ds["AR_INDEX"] = (
        ("time"),
        ar_events[1:-1].reshape(ds.IVT.shape),
    )
    ds["AR_DURATION"] = (
        ("time"),
        ar_duration[1:-1].reshape(ds.IVT.shape),
    )

    return ds


def compute_ar_categories(ds):
    """
        This function computes AR categories based on the AR index and AR duration.

    Parameters:
    ===========
        ds: xarray Dataset containing AR index and AR duration

    Returns:
    ========
        xarray dataset with AR categories
    """

    ds["AR_CATEGORY"] = (["time"], np.zeros(len(ds.AR_INDEX)))

    # compute AR Category by combining max_ivt and AR duration.
    for ar_id, d in ds.groupby(ds["AR_INDEX"]):

        max_ivt = int(round(d.IVT.max().values / 10, 0))
        duration = (
            int(d.AR_DURATION.max().values) * 6
        )  # don't need max here, but we do need a single value so this works for now.

        cat = AR[max_ivt, duration]

        ds["AR_CATEGORY"] = ds["AR_CATEGORY"].where(ds["AR_INDEX"] != ar_id, other=cat)
    return ds
