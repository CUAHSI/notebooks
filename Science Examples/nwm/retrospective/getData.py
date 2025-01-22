import os
import sys
import pytz
import urllib3
import datetime
import numpy as np
import pandas as pd
import pyproj
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
pd.options.mode.chained_assignment = None

def getData(SiteName, SiteID, StateAbb, StartDate, EndDate, OutputFolder):
	url1 = 'https://wcc.sc.egov.usda.gov/reportGenerator/view_csv/customMultiTimeSeriesGroupByStationReport/daily/start_of_period/'
	url2 = f'{SiteID}:{StateAbb}:SNTL%7Cid=%22%22%7Cname/'
	url3 = f'{StartDate},{EndDate}/'
	url4 = 'WTEQ::value?fitToScreen=false'
	url = url1+url2+url3+url4
	print(f'Start retrieving data for {SiteName}, {SiteID}')
	
	http = urllib3.PoolManager()
	response = http.request('GET', url)
	data = response.data.decode('utf-8')
	i=0
	for line in data.split("\n"):
		if line.startswith("#"):
			i=i+1
	data = data.split("\n")[i:]
	
	df = pd.DataFrame.from_dict(data)
	df = df[0].str.split(',', expand=True)
	df.rename(columns={0:df[0][0], 
					   1:df[1][0]}, inplace=True)
	df.drop(0, inplace=True)
	df.dropna(inplace=True)
	df.reset_index(inplace=True, drop=True)
	df["Date"] = pd.to_datetime(df["Date"])
	df.rename(columns={df.columns[1]:'Snow Water Equivalent (m) Start of Day Values'}, inplace=True)
	df.iloc[:, 1:] = df.iloc[:, 1:].apply(lambda x: pd.to_numeric(x) * 0.0254)  # convert in to m
	df['Water_Year'] = pd.to_datetime(df['Date']).map(lambda x: x.year+1 if x.month>9 else x.year)
	
	df.to_csv(f'./{OutputFolder}/df_{SiteID}_{StateAbb}_SNTL.csv', index=False)

def convert_latlon_to_yx(lat, lon, input_crs, ds, output_crs):
    """
    This function takes latitude and longitude values along with
    input and output coordinate reference system (CRS) and 
    uses Python's pyproj package to convert the provided values 
    (as single float values, not arrays) to the corresponding y and x 
    coordinates in the output CRS.
    
    Parameters:
    lat: The latitude value
    lon: The longitude value
    input_crs: The input coordinate reference system ('EPSG:4326')
    output_crs: The output coordinate reference system
    
    Returns:
    y, x: a tuple of the transformed coordinates in the specified output
    """
    # Create a transformer
    transformer = pyproj.Transformer.from_crs(input_crs, output_crs, always_xy=True)

    # Perform the transformation
    x, y = transformer.transform(lon, lat)

    return y, x 

def convert_utc_to_local(state_abbr, df):
    state_timezones = {
    'AL': 'US/Central', 'AK': 'US/Alaska', 'AZ': 'US/Mountain', 'AR': 'US/Central',
    'CA': 'US/Pacific', 'CO': 'US/Mountain', 'CT': 'US/Eastern', 'DE': 'US/Eastern',
    'FL': 'US/Eastern', 'GA': 'US/Eastern', 'HI': 'US/Hawaii', 'ID': 'US/Mountain',
    'IL': 'US/Central', 'IN': 'US/Eastern', 'IA': 'US/Central', 'KS': 'US/Central',
    'KY': 'US/Eastern', 'LA': 'US/Central', 'ME': 'US/Eastern', 'MD': 'US/Eastern',
    'MA': 'US/Eastern', 'MI': 'US/Eastern', 'MN': 'US/Central', 'MS': 'US/Central',
    'MO': 'US/Central', 'MT': 'US/Mountain', 'NE': 'US/Central', 'NV': 'US/Pacific',
    'NH': 'US/Eastern', 'NJ': 'US/Eastern', 'NM': 'US/Mountain', 'NY': 'US/Eastern',
    'NC': 'US/Eastern', 'ND': 'US/Central', 'OH': 'US/Eastern', 'OK': 'US/Central',
    'OR': 'US/Pacific', 'PA': 'US/Eastern', 'RI': 'US/Eastern', 'SC': 'US/Eastern',
    'SD': 'US/Central', 'TN': 'US/Central', 'TX': 'US/Central', 'UT': 'US/Mountain',
    'VT': 'US/Eastern', 'VA': 'US/Eastern', 'WA': 'US/Pacific', 'WV': 'US/Eastern',
    'WI': 'US/Central', 'WY': 'US/Mountain'
    }

    # Extract the state abbreviation from the filename
    # state_abbr = os.path.basename(filename).split('_')[2]  
    timezone = state_timezones.get(state_abbr)

    if timezone:
        # Convert the 'Date' column to datetime
        df['Date'] = pd.to_datetime(df['Date'], utc=True)
        
        # Convert to local time zone
        local_tz = pytz.timezone(timezone)
        df['Date_Local'] = df['Date'].dt.tz_convert(local_tz)

         # Save the timezone-aware Date_Local column
        df['Date_Local'] = df['Date_Local'].astype(str)
        df['Date_Local'] = df['Date_Local'].apply(lambda x: datetime.datetime.strptime(x, '%Y-%m-%d %H:%M:%S%z'))
        df['Date_Local'] = df['Date_Local'].apply(lambda x: x.replace(tzinfo=None))

    else:
        print(f"Timezone for state abbreviation {state_abbr} not found.")
        
    return df

def combine(snotel_files, nwm_files, StartDate, EndDate):

    # Create a dictionary to store dataframes
    dataframes = {}
    
    # Read SNOTEL files
    for file in snotel_files:
        location = os.path.basename(file).split('_')[1]  # Extract location from filename
        df = pd.read_csv(file)
        df['Date'] = pd.to_datetime(df['Date']).dt.date  # .dt.date is required because times are not excatly the same between SNOTEL and NWM
        dataframes[f'snotel_{location}'] = df.set_index('Date')
    
    # Read NWM files
    for file in nwm_files:
        location = os.path.basename(file).split('_')[1]  # Extract location from filename
        df = pd.read_csv(file)
        df['Date_Local'] = pd.to_datetime(df['Date_Local']).dt.date  # .dt.date is required because times are not excatly the same between SNOTEL and NWM
        dataframes[f'nwm_{location}'] = df.set_index('Date_Local')
    
    # Merge dataframes on Date
    combined_df = pd.DataFrame(index=pd.date_range(start=StartDate, end=EndDate))  
    for key, df in dataframes.items():
        if 'snotel' in key:
            combined_df[f'{key}_swe_m'] = df['Snow Water Equivalent (m) Start of Day Values']
        elif 'nwm' in key:
            combined_df[f'{key}_swe_m'] = df['NWM_SWE_meters']

    return combined_df

if __name__ == "__main__":
	SiteName = sys.argv[1]
	SiteID = sys.argv[2]
	StateAbb = sys.argv[3]
	StartDate = sys.argv[4]
	EndDate = sys.argv[5]
	OutputFolder = sys.argv[6]
	
	getData(SiteName, SiteID, StateAbb, StartDate, EndDate, OutputFolder)