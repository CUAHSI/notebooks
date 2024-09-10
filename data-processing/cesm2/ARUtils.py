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
            if flattened[idx - 1] == 0:
                # this is the first value in the set
                idx_start = idx
                event += 1
            current_sequence += 1
            ar_events[idx] = event

        # we hit a value of 0
        else:
            if current_sequence > 0:
                ar_duration[idx_start:idx] = current_sequence
            current_sequence = 0

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


######                                              #####
## Functions that use numpy arrays and are designed to ##
## be used with xarray.apply_ufunc                     ##
######                                              #####
def identify_ar_events(da):
    AR_THRESHOLD = 250

    # isolate the days in which the minimum AR conditions are met.
    ds_thr = xr.where((da >= AR_THRESHOLD), 1, 0)

    # Flatten the array
    flattened = ds_thr
    flattened = np.insert(flattened, 0, 0, axis=0)
    flattened = np.append(flattened, 0)

    ar_events = np.zeros(flattened.shape)

    # Initialize variables
    event = 0

    # Initialize variables
    for idx in range(1, len(flattened)):

        # we hit a value of 1
        if flattened[idx] == 1:
            if flattened[idx - 1] == 0:
                # this is the first value in the set
                event += 1
            ar_events[idx] = event

    # return all but the first zero that we added
    return ar_events[1:-1]


def compute_ar_durations(arr):
    """
    arr: DataArray of AR_INDEX
    """

    # Create an empty list to store the result
    result = []

    # Initialize variables to keep track of the current sequence
    current_value = arr[0]
    count = 1

    # Iterate over the array starting from the second element
    for i in range(1, len(arr)):

        if (
            arr[i] == arr[i - 1]
        ):  # Check if the current value is identical to the previous one
            # If the current number is identical, increment the count
            count += 1
        else:
            # skip values that equal zero
            if current_value == 0:
                result.extend([0] * count)
            else:
                # If the current number is not identical, append the count to the result list
                result.extend([count] * count)

            # Reset the current value and count
            current_value = arr[i]
            count = 1

    # Append the count for the last group of identical integers
    if current_value == 0:
        result.extend([0] * count)
    else:
        result.extend([count] * count)

    return np.array(result)


def compute_ar_category(ar_index, ar_duration, ar_ivt):

    # Initialize AR Category lookup table
    AR = np.zeros((200, 300)) + 5
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

    # initialize array to store category information
    cat = np.zeros(len(ar_index))

    # Initialize variables to keep track of the current sequence
    last_idx = 0
    max_ivt = 0
    duration = 0
    start_i = 0

    # compute AR Category by combining max_ivt and AR duration.
    # Iterate over the array starting from the second element
    for i in range(0, len(ar_index)):
        idx = ar_index[i]

        # condition 1: zero
        if idx == 0:
            if max_ivt > 0:
                cat[start_i:i] = AR[int(round(max_ivt / 10, 0)), int(duration * 6)]
                max_ivt = 0
                duration = 0

        # condition 2: current index is not zero and last index is zero
        elif idx > 0 and last_idx == 0:
            start_i = i
            duration = ar_duration[i]

        # condition 3: current index == last index
        elif idx == last_idx:
            # save the ivt if it's greater than max_ivt
            if ar_ivt[i] > max_ivt:
                max_ivt = ar_ivt[i]

        # save this index and move onto the next, i.e. ignore zeros.
        last_idx = idx

    return cat
