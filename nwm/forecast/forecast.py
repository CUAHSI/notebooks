#!/usr/bin/env python3

import time
import gcsfs
import shutil
import fsspec
import fnmatch
import xarray as xr
from tqdm import tqdm
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


def download_matching_files(bucket_name, date, forecast_mode, init_time, destination_folder, merge_files=False, clean_on_success=False, merge_format='Zarr'):
    """
    Download files from Google Cloud Storage that match a prefix.

    Args:
        bucket_name (str): The name of the GCS bucket.
        prefix (str): The prefix to match files against.
        wildcard (str): The wildcard pattern to filter files by.
        destination_folder (str): The local folder to save downloaded files.
        merge_files (bool): Indicates if a merged file will be created from the downloaded content.
        clean_on_success (bool) : Indicates if the downloaded files will be downloaded after merge is successgful.
        merge_format (str): Format to save the merged file 'Zarr' or 'NetCDF'.
    Returns:
        None
    """

    # Construct the prefix for the streamflow file
    prefix = f"nwm.{date.strftime('%Y%m%d')}/{forecast_mode}"
    Path(f'{destination_folder}/{prefix}').mkdir(parents=True, exist_ok=True)
    wildcard = f"nwm*t{init_time}z*channel*"

    # add the correct file extension to the merged_path
    merged_path = f'{destination_folder}/nwm.{date.strftime("%Y%m%d")}/t{init_time}z_{forecast_mode}'
    if merge_format == 'NetCDF':
        merged_path += '.nc'
        engine = 'h5netcdf'
    elif merge_format == 'Zarr':
        merged_path += '.zarr'
        engine = 'zarr'
    
    if Path(merged_path).exists():
        # exit early, no need to collect data
        print(f'  - Data already exists at {merged_path}, skipping download')
        return [merged_path]
        
    blobs = get_matching_blobs(bucket_name, f'{prefix}/', wildcard)
    print(f"  - Found {len(blobs)} matching files.")

    
    local_files = []
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = []

        # Download matching files
        for blob in blobs:
            # Construct local file path
            local_file_path = f"{destination_folder}/{blob.name}"
            
            local_files.append(local_file_path)
            
            # check if file already exists
            if Path(local_file_path).exists():
                #print(f"File {local_file_path} already exists. Skipping download.")
                continue
            future = executor.submit(download_blob, blob, local_file_path)
            futures.append(future)

        # progress bar for downloads
        tqdm_prefix = f'Downloading {forecast_mode} {date.strftime("%Y%m%d")}'
        with tqdm(total=len(futures), desc=tqdm_prefix) as pbar:
            for future in concurrent.futures.as_completed(futures):
                pbar.update(1)
    
    if merge_files:
        st = time.time()
        print('+ Merging files')
        
        ## TODO: CDO isn't working properly
        # from cdo import Cdo
        # cdo = Cdo()
        # cdo.cat(input=local_files, output=merged_path)

        # this is inefficient but works fine for now
        print('   - Reading NetCDF into Memory...', end='')
        ds = xr.open_mfdataset(local_files,
                               engine=engine,
                               parallel=True,
                               preprocess=lambda ds: ds[['time', 'streamflow', 'feature_id']])
        print('done')

        if merge_format == 'NetCDF':
            print('   - Saving to NetCDF...', end='')
            ds.to_netcdf(merged_path)
            print('done')
        elif merge_format == 'Zarr':
            print('   - Saving to Zarr...', end='')
            ds.to_zarr(merged_path)
            print('done')
        else:
            print(f'!! Unrecognized merge_format: {merge_format}, skipping')
        
        del ds
        print(f'  elapsed time {time.time() - st}')

        if clean_on_success:
            st = time.time()
            print('+ Cleaning files...', end='')
            
            path_to_delete = Path(local_file_path).parent
            shutil.rmtree(path_to_delete)
            print(f'  elapsed time {time.time() - st}')
            
        return [merged_path]

    return local_files


def get_streamflow_for_reaches(date, init_times=[], reach_ids=[], forecast_mode='medium_range_mem1', destination_folder=Path('.cache'), merge_files=False, clean_on_success=False, merge_format='Zarr'):
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

    dats = {}
    for it in init_times:
        
        # zero pad the initialization time
        init_time = str(it).zfill(2)
    
        # # Construct the prefix for the streamflow file
        # prefix = f"nwm.{date.strftime('%Y%m%d')}/{forecast_mode}"
        # wildcard = f"nwm*t{init_time}z*channel*"
    
        # create .cache directory if it doesn't already exist
        Path(".cache").mkdir(parents=True, exist_ok=True)
        
        # Download the streamflow files
        print('+ Collecting streamflow data...')
        paths = download_matching_files("national-water-model",
                                        date,
                                        forecast_mode,
                                        init_time,
                                        destination_folder,
                                        merge_files,
                                        clean_on_success,
                                        merge_format=merge_format)
    
        # # create list of paths that match .cache/prefix/wildcard
        # paths = list((Path(".cache")/prefix).glob(wildcard))
    
        if len(paths) == 0:
            raise ValueError(f"No files found matching prefix {prefix} and wildcard {wildcard}")
            return None

        label = f"{date.strftime('%Y%m%d')}-{forecast_mode}-t{it}z"
        # Load the streamflow data
        st = time.time()
        import pdb; pdb.set_trace() 
        if len(paths) == 1:
            print('+ Loading single-file streamflow data...', end='')
            ds = xr.open_dataset(paths[0]).sel(feature_id=reach_ids)
            dats[label] = ds
        else:
            print('+ Loading multi-file streamflow data...', end='')
            #streamflow_data = xr.open_mfdataset(paths)
            ds = xr.open_mfdataset(paths,
                                   engine='h5netcdf',
                                   parallel=True,
                                   preprocess=lambda ds: ds[['time', 'streamflow', 'feature_id']].sel(feature_id=reach_ids))    
            dats[label] = ds
        print(f'Elapsed time: {time.time()-st}')


    return dats

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
    st = time.time()
    ds = get_streamflow_for_reaches(datetime(2024, 4, 20), 0, 6076039, )
    print(f'Elapsed time: {time.time()-st}')

#    # open data directory from GCP without downloading.
#    from dask.distributed import Client
#    client = Client()
#    print(client.dashboard_link)
#    import time; st = time.time()
#    ds = open_mfdataset_from_gcp(bucket_name, prefix, wildcard, [6076039])
#    print(f'Elapsed time: {time.time()-st}')
#    import pdb; pdb.set_trace()






