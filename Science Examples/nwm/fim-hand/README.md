# Basic Flood Inundation Mapping using Height Above Nearest Drainage

**Authors:**  
   - Tony Castronova <acastronova@cuahsi.org>    
   - Irene Garousi-Nejad <igarousi@cuahsi.org>  
    
**Last Updated:** 05.22.2024

**Description**:  

The purpose of this Jupyter Notebook is to demonstrate the process of generating flood inundation maps (FIM) using pre-computed Height Above Nearest Drainage (HAND) raster maps. This process has been developed by the NOAA Office of Water Prediction, see the [inundation-mapping project](https://github.com/NOAA-OWP/inundation-mapping) for more information. 

There are two approaches to using the data provided by NOAA OWP to compute FIM; basic (simplified) mapping and mosaic mapping. The latter represents the state of practice in this domain, however for simplicity this notebook will demonstrate the former. After understanding the simplified approach it will be clear how to extend this work to the more complex *mosaic* methodology.

The FIM approach outlined in this notebook requires several input datasets. These can be obtained from the ESIP-hosted cloud store using the following commands (note: full instructions are provided in repository linked above):

```
aws s3 ls s3://noaa-nws-owp-fim/hand_fim/  

aws s3 sync s3://noaa-nws-owp-fim/hand_fim/outputs/fim_4_4_0_0/12090301 \
    /your_local_folder_name/12090301 
```

**Data Description**

- `hydroTable_0.csv`: pre-computed reach rating curves.
- `rem_zeroed_masked_0.tif`: pre-computed HAND raster grid.
- `gw_catchments_reaches_filtered_addedAttributes_crosswalked_0.gpkg`: added attributes for reaches.
- `demDerived_reaches_split_filtered_addedAttributes_crosswalked_0.gpkg`: DEM-derived reach geometries, used for visualization.

**Software Requirements**:  

The software and operating system versions used to develop this notebook are listed below. To avoid encountering issues related to version conflicts among Python packages, we recommend creating a new environment variable and installing the required packages specifically for this notebook.

Tested on: MacOS Ventura 13.2.1  

> boto3: 1.26.76  
  dask-core: 2023.4.0  
  fiona: 1.9.3  
  fsspec: 2023.4.0  
  geopandas: 0.12.2   
  ipyleaflet: 0.17.2  
  ipywidgets: 7.7.5   
  matplotlib: 3.7.1   
  netcdf4: 1.6.3   
  numpy: 1.24.2  
  pandas: 2.0.0  
  requests: 2.28.2  
  s3fs: 2023.4.0  
  scipy: 1.10.1  
  xarray: 2023.4.1
  
**Supplementary Code**

To simplify this notebook several *helper* functions have been develop that are referenced. These functions are located in a module called `nwm_utils`.