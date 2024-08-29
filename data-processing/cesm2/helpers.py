#!/usr/bin/env python3

import pandas
import numpy as np
import matplotlib.pyplot as plt

def compute_rolling_windows(dat, variable='RAIN', num_days = 45, ensembles=[0], center=False):
    # compute 45-day rolling windows for each of the ensemble members
    
    rolling_dict = {}
    for member_id in ensembles:
        print(f'Processing Member {member_id}...', end='')
        dat_mem = dat.isel(member_id=member_id)[variable]
        rolling_mem = (dat_mem.sum(dim=['lat', 'lon']).rolling(time=num_days, center=center).sum()).compute()
        rolling_dict[f'member-{member_id}'] = rolling_mem
        print('done') 
    return rolling_dict

def create_rolling_window_plots(rolling_dict,
                                title='Rolling Window',
                                ylabel='Variable',
                                xlabel=''):

    fig, ax = plt.subplots(3, 1, figsize=(10, 10))

    for member_id, rolling_mem in rolling_dict.items():
        rolling_mem.plot(ax=ax[0], label=member_id)
        rolling_mem.cumsum().plot(ax=ax[1], label=member_id)
    
    ax[0].set_title(f'{title}');
    ax[0].set_ylabel(ylabel);
    ax[0].legend(loc=(1.04, 0))
    
    ax[1].set_title(f'{title} - Cumulative');
    ax[1].set_ylabel(ylabel);
    ax[1].legend(loc=(1.04, 0))


    # convert DataArrays into a pandas DataFrame
    data = {}
    for k, v in rolling_dict.items():
        data[k] = v.values
    data['time'] = v.time
    df = pandas.DataFrame(data)
    
    # create box plot (using matplotlib directly so I can format the date axis)
    member_values = np.transpose(df.loc[:, df.columns != 'time'].to_numpy())
    time = df.time.to_numpy()
    
    ax[2].boxplot(member_values, boxprops={'color': 'lightblue'})
    ax[2].set_title(f'{title} - Box')
    ax[2].yaxis.grid(True)
    
    # create x ticks
    xticks = [tick for tick in range(0, member_values.shape[1], 20)]
    xlabels = [pandas.to_datetime(str(time[tick])).strftime('%m-%d-%Y') for tick in xticks]
    ax[2].set_xticks(xticks, xlabels)
    ax[2].set_ylabel(ylabel)
    
    plt.xticks(rotation=45)    
    plt.tight_layout()
    

def has_rain(dat, rolling_window_days=45, threshold=0):
    
    # lazy compute has_rain
    has_rain = (dat > threshold).astype(int)

    # compute rolling window
    rolling = has_rain.rolling(time=rolling_window_days).sum()

    # convert the datetime index into pandas datetime.
    rolling['time'] = rolling.indexes['time'].to_datetimeindex()

    return rolling
    