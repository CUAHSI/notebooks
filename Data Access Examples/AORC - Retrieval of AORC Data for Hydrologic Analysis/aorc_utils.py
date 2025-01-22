import xarray as xr
import fsspec
import pyproj

def get_conus_bucket_url(variable_code):
    """
    This function returns the S3 bucket url path for the CONSUS forcing data (in zarr format) for a given variable code.
    
    Parameters:
    variable_code (str): The code of the variable.

    Returns:
    str: S3 bucket url path for forcing data.
    """
    conus_bucket_url = f's3://noaa-nwm-retrospective-3-0-pds/CONUS/zarr/forcing/{variable_code}.zarr'
    return conus_bucket_url

def load_dataset(conus_bucket_url):
    """
    This function loads the dataset from the given S3 bucket url path.
    
    Parameters:
    conus_bucket_url (str): The S3 bucket URL path of the dataset to load.

    Returns:
    xr.Dataset: The loaded dataset as a xarray dataset
    """
    ds = xr.open_zarr(
        fsspec.get_mapper(
            conus_bucket_url,
            anon=True
        ),
        consolidated=True
    )
    return ds

def reproject_coordinates(ds, lon, lat, input_crs='EPSG:4326'):
    """
    This function reprojects the given lon and lat coordinates to the dataset's CRS.
    
    Parameters:
    ds (xr.Dataset): The xarray dataset containing the CRS information.
    lon (float): The longitude to reproject.
    lat (float): The latitude to reproject.
    input_crs (str): The CRS of the input coordinates.

    Returns:
    tuple: The reprojected coordinates (x, y).
    """
    output_crs = pyproj.CRS(ds.crs.esri_pe_string)  # CRS of AORC dataset (lcc)
    transformer = pyproj.Transformer.from_crs(input_crs, output_crs, always_xy=True)
    x, y = transformer.transform(lon, lat)
    return x, y

def get_aggregation_code(aggr_name):
    """
    Gets a aggregation code for a given aggregation name.

    Parameters:
    aggr_name (str): Name of the aggregation
    
    Returns:
    str: A code for aggregation
    """
    agg_options = {
        'hour':'h',
        'day':'d',
        'month':'ME',
        'year':'YE'
    }
    
    if aggr_name not in agg_options:
        raise Exception(f"{aggr_name} is not a valid aggregation name")
    
    return agg_options[aggr_name]

def get_variable_code(variable_name):
    """
    Gets a code for a given variable name for which data can be retrieved

    Parameters:
    variable_name (str): Name of the variable

    Returns:
    str: A variable code
    """
    variables = {
        'Total Precipitation':'precip', 
        'Air Temperature':'t2d', 
        'Specific Humidity':'q2d', 
        'Downward Long-Wave Radiation Flux':'lwdown',
        'Downward Short-Wave Radiation Flux':'swdown',
        'Pressure':'psfc',
        'U-Component of Wind':'u2d',
        'V-Component of Wind':'v2d'
    }

    if variable_name not in variables:
        raise Exception(f"{variable_name} is not a valid variable name")
    
    return variables[variable_name]

def get_time_code(interval_name):
    """
    Gets a time code for a given time interval name

    Parameters:
    interval_name (str): Name of the interval

    Returns:
    str: A time code
    
    """    
    time_attrs = {
        'hour':'h',
        'day':'d',
        'month':'M',
        'year':'Y'
    }
    if interval_name not in time_attrs:
        raise Exception(f"{interval_name} is not a valid interval name")
    
    return time_attrs[interval_name]