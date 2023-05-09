## Comparing basin averaged AORC precipitation datasets over the Logan River Watershed

### Content of this resource:
This example demonstrates how to collect and visualize the AORC (version 1.1) and PRISM precipitation data and calculate the basin averaged precipitation from both datasets over the Logan River Watershed in Utah (Figure below). This HydroShare resource contains the following items:

- **logan-watershed-huc12.shp**: This is the shapefile of the Logan River Watershed that contains seven features (HUC 12 catchments). This file is a subset of the National Hydrography Dataset (NHD) and was retrieved from the USGS website. 

- **Daily_PRISM_Precipitation_in_2010_for_LoganRiverWatershed.csv**: This is the PRISM daily precipitation spatially averaged across the Logan River Watershed. This file was created using *query-daily-prism-precipitation.ipynb*.  

- **query-daily-prism-precipitation.ipynb**: This is a computational Jupyter Notebook in Python that downloads the PRISM precipitation data and clipped them for the Logan River Watershed. 

- **query-aorc-and-compare-with-prism.ipynb**: This is a computational Jupyter Notebook in Python that retrieves the AORC forcing data, clipped them for the Logan River Watershed, and compare these data with the PRISM precipitation (*Daily_PRISM_Precipitation_in_2010_for_LoganRiverWatershed.csv*).

- **readme.md**: This is the readme file explaining the resource and its files.

- **case-study-logan-river-watershed.png**: This is an image showing the geographical location of the study site. This image is referenced in this readme file.

### Data used for this exercise:

*Analysis of Record for Calibration (AORC)*  
>The Analysis of Record for Calibration (AORC) is a gridded record of near-surface weather conditions covering the continental United States and Alaska and their hydrologically contributing areas. It is defined on a latitude/longitude spatial grid with a mesh length of ~800 m (30 arc seconds), and a temporal resolution of one hour. Elements include hourly total precipitation, temperature, specific humidity, terrain-level pressure, downward longwave and shortwave radiation, and west-east and south-north wind components. It spans the period from 1979 at Continental U.S. (CONUS) locations / 1981 in Alaska, to the near-present (at all locations). This suite of eight variables is sufficient to drive most land-surface and hydrologic models and is used to force the calibration run of the National Water Model (NWM).

*Parameter-elevation Regressions on Independent Slopes Model (PRISM)*
> The PRISM Climate Data from Oregon State University is a gridded product covering the Conterminous United States. It combines climate observations from a wide range of monitoring networks with sophisticated modeling techniques to establish a local climate-elevation relationship for each grid cell. This relationship is then utilized to estimate various climate variables such as precipitation. This dataset spans from 1895 to the present and can be used to analyze short- and long-term climate patterns. The PRISM dataset is widely used across multiple disciplines, such as hydrology, climatology, agriculture, and environmental sciences and its popularity stems from its high spatial resolution, long-term record, reliable data quality, and ease of accessibility. It is worth noting that the PRISM dataset is available at two distinct spatial scales. The 4 km version is accessible to the public, while a more refined 800 m resolution dataset is available for a fee. For our purposes, we will be utilizing the publicly available 4 km resolution dataset.

*Watershed Boundary Dataset (WBD)*:
> The Watershed Boundary Dataset (WBD) is a seamless, national hydrologic unit dataset. Hydrologic units represent the area of the landscape that drains to a portion of the stream network. More specifically, a hydrologic unit defines the areal extent of surface water drainage to an outlet point on a dendritic stream network or to multiple outlet points where the stream network is not dendritic. WBD contains eight levels of progressive hydrologic units (HUC) identified by unique 2- to 16-digit codes. The case study region was selected from the HUC 12 watershed boundaries for the Great Basin from [this](https://www.hydroshare.org/resource/965eab1801c342a58a463f386c9f3e9b/) HydroShare resource. 

### Case study:
The Logan River Watershed is a snowmelt-dominated watershed located in the Bear River mountain range east of Logan, Utah. Most precipitation falls in the form of snowfall during winter months and as rain during spring and summer time. Water flows southwest through mostly natural land cover. 

<img src="https://www.hydroshare.org/resource/709bc880194046e3b79b0da56ad090fa/data/contents/case-study-logan-river-watershed.png" alt="Logan River Watershed" width="750" height="300">


## How to run computational notebooks:
To access the CUAHSI JupyterHub, you have two options. First, you can right-click on any of the Jupyter Notebooks within this resource and choose the `CUAHSI JupyterHub` option. Alternatively, you can open the entire resource by clicking on the `Open with` button located at the top right corner of the landing page. Once you're in the resource, make sure to select the `Python 3.8` server, and then follow the steps provided in the notebook.
