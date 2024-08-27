#!/usr/bin/env python3


def has_rain(dat, rolling_window_days=45, threshold=0):
    
    # lazy compute has_rain
    has_rain = (dat > threshold).astype(int)

    # compute rolling window
    rolling = has_rain.rolling(time=rolling_window_days).sum()

    # convert the datetime index into pandas datetime.
    rolling['time'] = rolling.indexes['time'].to_datetimeindex()

    return rolling
    