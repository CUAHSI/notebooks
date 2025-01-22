#!/usr/bin.env python3

import gcsfs
import fsspec
import fnmatch
import xarray as xr
from pathlib import Path
import concurrent.futures
from datetime import datetime
from google.cloud import storage


def get_matching_blobs(bucket_name, prefix, wildcard):
    """
    Gets a list of blobs from Google Cloud Storage that match a prefix.

    Args:
        bucket_name (str): The name of the GCS bucket.
        prefix (str): The prefix to match files against.
        wildcard (str): The wildcard pattern to filter files by.
    Returns:
        List of matching blobs
    """

    # Create a client for interacting with GCS
    storage_client = storage.Client.create_anonymous_client()

    # Get the specified bucket
    bucket = storage_client.bucket(bucket_name)

    # List blobs (files) with the given prefix
    blobs = bucket.list_blobs(prefix=prefix)

    # filter filenames by wildcard
    blobs = [blob for blob in blobs if fnmatch.fnmatchcase(blob.name, wildcard)]

    return blobs


def download_blob(blob, local_file_path):
    blob.download_to_filename(local_file_path)
    return blob.name


def download_matching_files(bucket_name, prefix, wildcard, destination_folder):
    """
    Download files from Google Cloud Storage that match a prefix.

    Args:
        bucket_name (str): The name of the GCS bucket.
        prefix (str): The prefix to match files against.
        wildcard (str): The wildcard pattern to filter files by.
        destination_folder (str): The local folder to save downloaded files.
    Returns:
        None
    """

    blobs = get_matching_blobs(bucket_name, prefix, wildcard)
    print(f"Found {len(blobs)} matching files.")

    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = []

        # Download matching files
        for blob in blobs:
            # Construct local file path
            local_file_path = f"{destination_folder}/{blob.name}"
            
            # check if file already exists
            if Path(local_file_path).exists():
                #print(f"File {local_file_path} already exists. Skipping download.")
                continue
            future = executor.submit(download_blob, blob, local_file_path)
            futures.append(future)

        for future in concurrent.futures.as_completed(futures):
            print(f"Downloaded {future.result()}")


def get_streamflow_for_reaches(date, init_time, reach_ids=[], forecast_mode='medium_range_mem1', destination_folder=Path('.cache')):
    """
    Returns a pandas dataframe of streamflow for a given reach and date.

    Parameters
    ----------
    date : datetime
        The date for which to retrieve streamflow data.
    init_time: int
        The initialization time for the streamflow data.
    reach_id : int
        The reach ID for which to retrieve streamflow data.
    """

    # zero pad the initialization time
    init_time = str(init_time).zfill(2)

    # Construct the prefix for the streamflow file
    prefix = f"nwm.{date.strftime('%Y%m%d')}/{forecast_mode}/"
    wildcard = f"nwm*t{init_time}z*channel*"


    # create .cache directory if it doesn't already exist
    Path(".cache").mkdir(parents=True, exist_ok=True)
    
    # Download the streamflow files
    print('+ Collecting streamflow data...')
    download_matching_files("national-water-model", prefix, wildcard, destination_folder)

    # create list of paths that match .cache/prefix/wildcard
    paths = list((Path(".cache")/prefix).glob(wildcard))

    if len(paths) == 0:
        raise ValueError(f"No files found matching prefix {prefix} and wildcard {wildcard}")
        return None

    # Load the streamflow data
    print('+ Loading streamflow data...', end='')
    st = time.time()
    streamflow_data = xr.open_mfdataset(paths)
    print(f'Elapsed time: {time.time()-st}')

    # Extract the streamflow for the given reach
    streamflow = streamflow_data.streamflow.sel(feature_id=reach_ids)

    return streamflow

# TODO: this is super slow. Need to figure out how to speed up the download process.
def open_mfdataset_from_gcp(bucket_name, prefix, wildcard, reach_ids):
    """
    Open multiple files as a single dataset.
    Args:
        bucket_name (str): The name of the GCS bucket.
        prefix (str): The prefix to match files against.
        wildcard (str): The wildcard pattern to filter files by.
        reach_ids (List(int)): List of NWM reach ids to return data for
    Returns:
        xarray.Dataset
    """
    
    # connect to gcs anonymously
    fs = gcsfs.GCSFileSystem(token=None)

    # list files that match the prefix and wildcard
    files = fs.glob(f'{bucket_name}/{prefix}/{wildcard}')
    
    # open the files as a single dataset
    ds = xr.open_mfdataset([fs.open(f) for f in files[:25]],
                           engine='h5netcdf',
                           parallel=True,
                           preprocess=lambda ds: ds[['time', 'streamflow', 'feature_id']].sel(feature_id=reach_ids))


    # return the dataset with only the specified reach ids
    return ds 
    

# Example usage
if __name__ == "__main__":
    bucket_name = "national-water-model"
    prefix = "nwm.20240420/medium_range_mem1"
    wildcard = "nwm.t00z*channel_rt*"
    destination_folder = ".cache"
    

    # create output path if it doesn't already exist
    output_path = Path(destination_folder)/prefix
    output_path.mkdir(parents=True, exist_ok=True)

#    wildcard = "nwm*t00z*channel_rt*" # periods cannot be included for fnmatch to work
#    download_matching_files(bucket_name, prefix, wildcard, destination_folder)

    # collect data for the Squannacook River at West Groton
    import time; st = time.time()
    ds = get_streamflow_for_reaches(datetime(2024, 4, 20), 0, 6076039)
    print(f'Elapsed time: {time.time()-st}')

#    # open data directory from GCP without downloading.
#    from dask.distributed import Client
#    client = Client()
#    print(client.dashboard_link)
#    import time; st = time.time()
#    ds = open_mfdataset_from_gcp(bucket_name, prefix, wildcard, [6076039])
#    print(f'Elapsed time: {time.time()-st}')
#    import pdb; pdb.set_trace()






