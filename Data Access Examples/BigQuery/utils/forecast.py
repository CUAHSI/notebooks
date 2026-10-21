#!/usr/bin/env python3

"""
Description: This script contains helper functions for collecting National 
             Water Model timeseries data using the BiGQuery API, located 
             at https://nwm-api.ciroh.org/

Author(s): Tony Castronova <acastronova@cuahsi.org>
"""

import io
import pandas
import requests
import concurrent.futures
import matplotlib.pyplot as plt
from concurrent.futures import ThreadPoolExecutor


from utils import nwmmap
from datetime import datetime, timedelta
from ipyleaflet import Map, basemaps, WidgetControl
from ipywidgets import IntSlider, ColorPicker, jslink, Label, Text
import ipywidgets as widgets
from datetime import datetime, timedelta


# use the notebook version of TQDM if running in Jupyter
# otherwise use the regular version.
try:
    from tqdm.notebook import tqdm as tqdm
except ImportError:
    from tqdm import tqdm

                
class Forecasts():
    def __init__(self, api_key, api_url='https://nwm-api.ciroh.org'):
        self.url = f'{api_url}/forecast'
        self.header =  {
            'x-api-key': api_key
        }
        self.df = None
        
    
    def fetch_url(self, params):
        # try:
        response = requests.get(self.url,
                                params=params,
                                headers=self.header,
                                timeout=60)
        
        # Raise an exception for HTTP errors
        response.raise_for_status()  
        return response
            
        # except requests.exceptions.RequestException as e:
        #     return f"Error fetching {self.url}: {e}"

    def fetch_async(self, params_list):

        results = []
        errors = []
        
        # Use ThreadPoolExecutor to make concurrent GET requests
        # TQDM is used to provide a nice looking progress bar
        with ThreadPoolExecutor(max_workers=5) as executor:
            
            # Submit all URLs to the executor
            future_to_url = {executor.submit(self.fetch_url, param): param for param in params_list}
            
            # Process the results as they complete
            for future in tqdm(concurrent.futures.as_completed(future_to_url),
                               total=len(future_to_url),
                               desc="Fetching Forecast Data",
                               unit="url",
                               colour="green",  
                               dynamic_ncols=True):
                url = future_to_url[future]
                try:
                    res = future.result()

                    # If the request wasn't successful,
                    # log it as an error.
                    res.raise_for_status()

                    # otherwise, the 
                    results.append(res)
                    
                except Exception as e:
                    errors.append(
                        {
                            'error_message': e,
                            'parameters': params_list
                        }
                    )
            
            return results, errors
            
        
            
    def collect_forecasts(self,
                          comids, 
                          forecast_type,
                          ensemble,
                          reference_times):

        # build a parameters to query
        params = [
            {'comids': ','.join(map(str, comids)), 
             'forecast_type': forecast_type,
             'reference_time': reftime,
             'ensemble': ','.join(map(str, ensemble)),
             'output_format': 'csv'}
            for reftime in reference_times
        ]

        # query the api asynchronously with the parameters defined above 
        responses, errors = self.fetch_async(params)

        if errors:
            print(f'Encountered {len(errors)} errors')
            
        # filter out only the successful responses and 
        # convert them into a single pandas dataframe
        successful_responses = [resp for resp in responses if resp.status_code == 200]
        dfs = [pandas.read_csv(io.StringIO(res.text), sep=',') for res in successful_responses]
        df = pandas.concat(dfs, ignore_index=True)  

        # clean datetime columns and return
        df.time = pandas.to_datetime(df.time)
        df.reference_time = pandas.to_datetime(df.reference_time)

        self.df = df

    def plot(self, comid, plot_type='series', xlabel='Time', ylabel='Streamflow'):
        
        if self.df is None:
            print('No forecast data to plot. Run "collect_forecasts" to collect data')
            return None
            
        df = self.df[self.df.feature_id == int(comid)]
        
        fig, ax = plt.subplots(figsize=(10, 5))
        
        if plot_type == 'series':
            self.__plot_series(df, ax)
        elif plot_type == 'iqr':
            self.__plot_iqr(df, ax)
        else:
            print(f'Unrecognized plot type: {plot_type}')
            return None

        ax.grid(True)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        plt.xticks(rotation=25) 
        
        return plt, ax

    def __plot_series(self, df, ax):
        
        
        # Group by 'reference_time' and plot each group on the same axis
        for reference_time, group in df.groupby('reference_time'):
            group.plot(x='time', y='streamflow', ax=ax, label=str(reference_time), legend=False)
            
    def __plot_iqr(self, df, ax):
        
        iqr = df.groupby(df.time)['streamflow'].quantile([0.25, 0.75])
        iqr = iqr.reset_index()
        iqr = iqr.rename(columns={'level_1': 'quantile'})
        
        df_pivot = iqr.pivot(index='time', columns='quantile', values='streamflow')
        df_pivot.index = pandas.to_datetime(df_pivot.index)

        ax.fill_between(df_pivot.index, df_pivot[0.25], df_pivot[0.75], color='blue', alpha=0.3);


class ForecastMap(nwmmap.Map):
    def __init__(self, api_key):
        # initialize the parent class
        super().__init__()
        
        self.reach_label = None
        self.forecasts = Forecasts(api_key)

        
        # Plot Widget
        self.output = widgets.Output()
        self.output_close = widgets.Button(description="Close")
        self.output_close.on_click(self.on_plot_close)
        self.plot_container = widgets.VBox([self.output, self.output_close])      
        self.plot_container.layout.display = 'none'
        self.plot_control = WidgetControl(widget=self.plot_container, position='bottomleft')
        self.map.add(self.plot_control)
        
        
        self.styled_container = widgets.HTML("<style>.floating-widget { z-index: 9999 !important; position: relative; }</style>")
        self.map.add(WidgetControl(widget=self.styled_container, position="bottomleft"))



        # Forecast Submit Widget        
        self.submit = widgets.Button(description='Submit',
                                     disabled=True,
                                     button_style='',
                                     icon=' ')
        self.submit.on_click(self.collect_forecasts)
        
        self.map.add(WidgetControl(widget=self.submit,
                                   position='bottomleft'))

        # Forecast Option Widgets
        self.reach_label = widgets.Text(value='No Reach Selected',
                                        disabled=True,  
                                        layout=widgets.Layout(width='200px'))
        
        self.start_date = widgets.DatePicker(disabled = False,
                                             value = datetime.today() - timedelta(days=10),
                                             layout=widgets.Layout(width='200px'))

        self.forecast_type = widgets.Dropdown(options=['short_range', 'medium_range', 'long_range'],
                                              value='medium_range',
                                              disabled=False,
                                              layout=widgets.Layout(width='200px'))
        
        self.num_days = widgets.IntText(value=10,
                                        disabled=False,
                                        layout=widgets.Layout(width='200px'))

        label_layout = widgets.Layout(width='100px', justify_content='flex-end')

        self.options_widgets = widgets.VBox([
            widgets.HBox([Label(value='Reach', layout=label_layout), self.reach_label]),
            widgets.HBox([Label(value='Start Date', layout=label_layout), self.start_date]),
            widgets.HBox([Label(value='Forecast Type', layout=label_layout), self.forecast_type]),
            widgets.HBox([Label(value='Number of Days', layout=label_layout), self.num_days])
        ])
        
        self.map.add(WidgetControl(widget=self.options_widgets,
                                   position='bottomleft'))


        # debugging - REMOVE ME
        self.reach_label.value = '4965151'
        self.submit.disabled = False
        

        
    def on_plot_close(self, callback):
        print('on_plot_close')
        self.output.clear_output(wait=True)
        
        # toggle the visibility of map widgets
        self.toggle_widget_visibility([self.plot_container, self.options_widgets, self.submit])
        

    def action_after_map_click(self):

        if self.selected() is not None:
            # set the value of the label
            self.reach_label.value = self.selected()

            # activate the submission button
            self.submit.disabled = False
            self.submit.button_style = 'info'
            
        else:
            # set the value of the label
            self.reach_label.value = 'No Reach Selected'

            # deactivate the submission button
            self.submit.disabled = True
            self.submit.button_style = ''

    def collect_forecasts(self, callback):
            
        reference_times = [(self.start_date.value + timedelta(days=i)).strftime('%Y-%m-%d')
                           for i in range(self.num_days.value)]
        
        self.forecasts.collect_forecasts([self.reach_label.value],
                                         self.forecast_type.value,
                                         [0],
                                         reference_times)
        self.show_plot()

    def show_plot(self):

        # toggle the visibility of map widgets
        self.toggle_widget_visibility([self.plot_container, self.options_widgets, self.submit])
        
        with self.output:            

            # create and show a plot of the data
            plt, ax = self.forecasts.plot(self.reach_label.value, plot_type='iqr', ylabel='Streamflow (cms)')
            ax.set_title(f'IQR of Forecasted Streamflow for NWM {self.reach_label.value}')
            plt.show()

    def toggle_widget_visibility(self, widgets=[]):

        for widget in widgets:
            visibility = widget.layout.display
            
            if visibility == 'none':
                widget.layout.display = 'flex'
            else:
                widget.layout.display = 'none'
            
            

