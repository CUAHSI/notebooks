# Authors: Tony Castronova acastronova@cuahsi.org, Irene Garousi-Nejad igarousi@cuahsi.org
# Last Updated: 04.04.2023
# This script contains two function that are based on the workflow previously published on https://github.com/CUAHSI/notebook-examples/

import os
import re
import numpy
import xarray
import pyproj
import requests
import rioxarray
from pyproj import Transformer
import xml.etree.ElementTree as ET
from owslib.wms import WebMapService


def get_paths(subset_times):
    '''
    Query the AORC v1.0 data on HydroShare's Thredds Server
    A Jupyter notebook verions of this function can be found at https://github.com/CUAHSI/notebook-examples/blob/main/thredds/query-aorc-thredds.ipynb
    Input: A list of monthly datetime strings that are within the period of interest.
    Output: A list of path to AORC v1.0 data on HydroShed's Thredds
    '''
    
    catalog_base_url = 'https://thredds.hydroshare.org/thredds/catalog'
    dods_base_url = 'https://thredds.hydroshare.org/thredds/dodsC'

    url = f'{catalog_base_url}/aorc/data/16/catalog.xml'
    root = ET.fromstring(requests.get(url).text)
    ns = '{http://www.unidata.ucar.edu/namespaces/thredds/InvCatalog/v1.0}'

    # use xpath top select all "dataset" elements.
    elems = root.findall(f'.//{ns}dataset')

    # loop through results and extract the "urlPath" attribute values
    paths = []
    for elem in elems:
        atts = elem.attrib
        if 'urlPath' in atts.keys():
            paths.append(f"{dods_base_url}/{atts['urlPath']}")

    # use regex to isolate only files that end with ".nc"
    paths = list(filter(re.compile("^.*\.nc$").match, paths))
    
    
    files = []
    for i in range(0, len(subset_times)):
        f = str(subset_times[i]).split(" ")[0].split("-")[0:2][0]+str(subset_times[i]).split(" ")[0].split("-")[0:2][1]+'.nc'
        files.append(f)

    # Create a set of file names for faster lookup
    file_names = set(files)
    
    # Create a list to store the paths that include file names
    selected_paths = []
    
    # Iterate over each path in paths
    for path in paths:
        # Extract the file name from the path
        file_name = path.split("/")[-1]
        # Check if the file name is in the set of file names
        if file_name in file_names:
            # Append the path to the selected_paths list
            selected_paths.append(path)
    
    return selected_paths 



def add_spatial_metadata(ds):
    '''
    Add GeoSpatial Metadata to ACOR v1.0 Forcing Data
    A Jupyter notebook verions of this function can be found at https://github.com/CUAHSI/notebook-examples/blob/main/thredds/aorc-adding-spatial-metadata.ipynb
    Input:
        ds: A dataset for which we want to add geospatial metadata. 
    Output: A dataset that includes both projected and geographical coordinates.
    '''
    
    # load the meta datafile
    ds_meta = xarray.open_dataset('http://thredds.hydroshare.org/thredds/dodsC/hydroshare/resources/2a8a3566e1c84b8eb3871f30841a3855/data/contents/WRF_Hydro_NWM_geospatial_data_template_land_GIS.nc')

    
    def pattern_lookup(pattern, input): 
    
        # use the re.search() function to search for the pattern in the string
        match = re.search(pattern, input)

        # check if a match was found
        if match:
            # extract the matched values and concatenate them into the desired string format
            result = f'{match.group(0)}'
            return result
        else:
            # if no match was found, print an error message
            print('No match found.')
        
        
    pattern_we = r'west_east,(\d+),(\d+)'
    pattern_sn = r'south_north,(\d+),(\d+)'
    
    GSL_westeast = pattern_lookup(pattern_we, ds.attrs['history'])
    GSL_southnorth = pattern_lookup(pattern_sn, ds.attrs['history'])

    y_index = GSL_southnorth.split(',')[1:]
    x_index = GSL_westeast.split(',')[1:]
    
    leny = len(ds_meta.y)
    x = ds_meta.x[int(x_index[0]) : int(x_index[1]) + 1].values
    y = ds_meta.y[leny - int(y_index[1]) - 1 : leny - int(y_index[0])].values
    
    ds = ds.rename_dims(south_north='y', west_east='x', Time='time')
    
    X, Y = numpy.meshgrid(x, y)

    # define the input crs
    wrf_proj = pyproj.Proj(proj='lcc',
                           lat_1=30.,
                           lat_2=60., 
                           lat_0=40.0000076293945, lon_0=-97., # Center point
                           a=6370000, b=6370000) 

    # define the output crs
    wgs_proj = pyproj.Proj(proj='latlong', datum='WGS84')

    # transform X, Y into Lat, Lon
    transformer = Transformer.from_crs(wrf_proj.crs, wgs_proj.crs)
    lon, lat = transformer.transform(X, Y)
    
    
    ds = ds.assign_coords(lon = (['y', 'x'], lon))
    ds = ds.assign_coords(lat = (['y', 'x'], lat))
    ds = ds.assign_coords(x = x)
    ds = ds.assign_coords(y = y)
    
    
    ds.x.attrs['axis'] = 'X'
    ds.x.attrs['standard_name'] = 'projection_x_coordinate'
    ds.x.attrs['long_name'] = 'x-coordinate in projected coordinate system'
    ds.x.attrs['resolution'] = 1000.  # cell size
    ds.x.attrs['units'] = 'm'

    ds.y.attrs['axis'] = 'Y' 
    ds.y.attrs['standard_name'] = 'projection_y_coordinate'
    ds.y.attrs['long_name'] = 'y-coordinate in projected coordinate system'
    ds.y.attrs['resolution'] = 1000.  # cell size
    ds.y.attrs['units'] = 'm'

    ds.lon.attrs['units'] = 'degrees_east'
    ds.lon.attrs['standard_name'] = 'longitude' 
    ds.lon.attrs['long_name'] = 'longitude'

    ds.lat.attrs['units'] = 'degrees_north'
    ds.lat.attrs['standard_name'] = 'latitude' 
    ds.lat.attrs['long_name'] = 'latitude'
    
    # add crs to netcdf file
    ds.rio.write_crs(ds_meta.crs.attrs['spatial_ref'], inplace=True
                    ).rio.set_spatial_dims(x_dim="x",
                                           y_dim="y",
                                           inplace=True,
                                           ).rio.write_coordinate_system(inplace=True)

    return ds    
        