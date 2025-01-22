
import re
import dask
import json
import numpy
import xarray
import pyproj
import pandas
import requests
import geopandas
from matplotlib import colors
import matplotlib.pyplot as plt
from dask.distributed import Client
from dask.distributed import progress

import zarr
import fsspec
from pyproj import Transformer
from s3fs import S3FileSystem
from kerchunk.combine import MultiZarrToZarr

import rioxarray
from geocube.api.core import make_geocube


# |%%--%%| <tf48s3BaYX|uviRDjtlSm>

# use a try accept loop so we only instantiate the client
# if it doesn't already exist.
try:
    print(client.dashboard_link)
except:    
    # The client should be customized to your workstation resources.
    client = Client(n_workers=2, memory_limit='2GB') # per worker
    print(client.dashboard_link)

# |%%--%%| <uviRDjtlSm|xJLgZw2nBH>
r"""°°°
## Load Forcing Data into Memory


In this notebook we'll be working with AORC v1.0 meteorological forcing. These data are publicly available for the entire CONUS, spanning from 1980 to 2020. Kerchunk header files have been created by the Alabama Water Institute team and this is an ongoing project. Please note that this jupyter notebook works for data within 2007-2019, but it cannot work with data prior to 2006.
°°°"""
# |%%--%%| <xJLgZw2nBH|FFKy9uZGuD>

# define the selected watershed boundary 
wb_id = 'wb-2851655'

# define the year of interest
year=2010

# |%%--%%| <FFKy9uZGuD|iLCRKslN8J>

bucket = 's3://ciroh-nwm-zarr-retrospective-data-copy/noaa-nwm-retrospective-2-1-zarr-pds/forcing/'

# create an instace of the S3FileSystem class from s3fs
s3 = S3FileSystem(anon=True)
files = s3.ls(f'{bucket}{year}')  

new_files = []
for f in files:
    parts = f.split('/')
    parts[0] += '.s3.amazonaws.com'
    parts.insert(0, 'https:/')
    new_name = '/'.join(parts)
    new_files.append(new_name)
    

# |%%--%%| <iLCRKslN8J|n1L354D1TE>

%%time
json_list = new_files[0:217] 

mzz = MultiZarrToZarr(json_list,
    remote_protocol='s3',
    remote_options={'anon':True},
    concat_dims=['valid_time'])

d = mzz.translate()

backend_args = {"consolidated": False, "storage_options": {"fo": d}, "consolidated": False}

ds = xarray.open_dataset("reference://", engine="zarr", backend_kwargs=backend_args)

ds = ds.squeeze(dim='Time')

# |%%--%%| <n1L354D1TE|Jez1yBTVcA>

# !wget https://thredds.hydroshare.org/thredds/fileServer/hydroshare/resources/2a8a3566e1c84b8eb3871f30841a3855/data/contents/WRF_Hydro_NWM_geospatial_data_template_land_GIS.nc

# |%%--%%| <Jez1yBTVcA|tosYKrYoiv>


ds_meta = xarray.open_dataset('WRF_Hydro_NWM_geospatial_data_template_land_GIS.nc')
leny = len(ds_meta.y)
x = ds_meta.x.values
y = ds_meta.y.values

ds = ds.rename({'valid_time': 'time', 'south_north':'y', 'west_east':'x'})

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
transformer = pyproj.Transformer.from_crs(wrf_proj.crs, wgs_proj.crs)
lon, lat = transformer.transform(X, Y)

ds = ds.assign_coords(lon = (['y', 'x'], lon))
ds = ds.assign_coords(lat = (['y', 'x'], lat))
ds = ds.assign_coords(x = x)
ds = ds.assign_coords(y = y)

ds.x.attrs['axis'] = 'X'
ds.x.attrs['standard_name'] = 'projection_x_coordinate'
ds.x.attrs['long_name'] = 'x-coordinate in projected coordinate system'
ds.x.attrs['resolution'] = 1000.  # cell size

ds.y.attrs['axis'] = 'Y' 
ds.y.attrs['standard_name'] = 'projection_y_coordinate'
ds.y.attrs['long_name'] = 'y-coordinate in projected coordinate system'
ds.y.attrs['resolution'] = 1000.  # cell size

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
                                       ).rio.write_coordinate_system(inplace=True);

# |%%--%%| <tosYKrYoiv|MCvi4hcKRC>
r"""°°°
## Add spatial reference to the model domain
°°°"""
# |%%--%%| <MCvi4hcKRC|EqPKq2ucbN>

# prepare geometries for spatial averaging
gdf = geopandas.read_file(f'domain/v20.1/{wb_id}/{wb_id.split("_")[0]}_upstream_subset.gpkg', layer='divides')


# convert these data into the projection of our forcing data
target_crs = pyproj.Proj(proj='lcc',
                       lat_1=30.,
                       lat_2=60., 
                       lat_0=40.0000076293945, lon_0=-97., # Center point
                       a=6370000, b=6370000) 

gdf = gdf.to_crs(target_crs.crs)

# important step
# rechunk the dataset to solve the memory limit issue
ds = ds.chunk(chunks={'time': 1})

# |%%--%%| <EqPKq2ucbN|CqHhLNEnyW>
r"""°°°
## Build Weights Variable
°°°"""
# |%%--%%| <CqHhLNEnyW|zAcXQDBJbI>

#!wget 'https://lynker-spatial.s3.amazonaws.com/v20.1/forcing_grids.json'

# |%%--%%| <zAcXQDBJbI|64C8zs8pqA>

def assert_grids_equal(ds, weights_grid):
    dsx = ds.x.resolution
    dsy = ds.y.resolution
    assert dsx == weights_grid[0]['resX']
    assert dsy == weights_grid[0]['resY']
    
    assert ds.x.min().item() - dsx / 2 == weights_grid[0]['X1']
    assert ds.x.max().item() + dsx / 2 == weights_grid[0]['Xn']
    assert ds.y.min().item() - dsy / 2 == weights_grid[0]['Y1']
    assert ds.y.max().item() + dsy / 2 == weights_grid[0]['Yn']
    
    
    ds_nrows, ds_ncols = ds.LWDOWN.shape[1:]
    assert ds_ncols   == weights_grid[0]['ncols']
    assert ds_nrows   == weights_grid[0]['nrows']

    #ds_wkt = ds.spatial_ref.crs_wkt
    #crs = wgrid[0]['crs']

    return True


# |%%--%%| <64C8zs8pqA|0hLJyYjY8T>

wgrid = json.load(open('forcing_grids.json', 'r'))
assert_grids_equal(ds, wgrid)

# |%%--%%| <0hLJyYjY8T|cVWggJdhEC>

# Since our data are the same shape, compute cell indice for all x,y coords

xs = ds.x.values
ys = ds.y.values
dx = ds.x.resolution
dy = ds.y.resolution
xmin = ds.x.min().item() - dx / 2 
xmax = ds.x.max().item() + dx / 2 
ymin = ds.y.min().item() - dy / 2 
ymax = ds.y.max().item() + dy / 2 
ds_nrows, ds_ncols = ds.LWDOWN.shape[1:]


# |%%--%%| <cVWggJdhEC|pkNncKmglw>

# unravel x, y coordinates into point pairs

rr, cc = numpy.meshgrid(xs, ys)
pair_pts = numpy.array([rr, cc]).T.reshape(-1, 2)

# |%%--%%| <pkNncKmglw|HBaNBhAl2y>

xss = pair_pts[:,0]
yss = pair_pts[:,1]

# |%%--%%| <HBaNBhAl2y|oNrLaPCaI8>

cell_number = ((ymax - (yss - dy/2)) / dy * ds_ncols + ((xss- dy/2) - xmin) / dx).round()

# |%%--%%| <oNrLaPCaI8|yNuI3LuLLX>

cell_ids = cell_number.reshape(ds_nrows, ds_ncols)
cell_ids.shape

# |%%--%%| <yNuI3LuLLX|uDqdutBNlX>

# add as a new variable to the dataset
ds['cell_number']=(['y','x'], cell_ids)

# |%%--%%| <uDqdutBNlX|9oykqAY2RW>

ds

# |%%--%%| <9oykqAY2RW|osnjikLInq>
r"""°°°
## Add catchment ids based on the geometries in the subsetted hydrofabric
°°°"""
# |%%--%%| <osnjikLInq|2t1NIZnt9Y>
r"""°°°
Isolate the data corresponding to the subsetted area.
°°°"""
# |%%--%%| <2t1NIZnt9Y|npVWRBHKNy>

%%time

# create zonal id column
gdf['cat'] = gdf.id.str.split('-').str[-1].astype(int)

# clip AORC to the extent of the hydrofabric geometries
# todo : buffer gdf.geometry.values to include neighboring cells
buffered = gdf.dissolve().buffer(1000)
ds_subset = ds.rio.clip(buffered.geometry.values,
                 gdf.crs,
                 drop=True,
                 invert=False,
                 from_disk=True,
                 all_touched=True
                )

# |%%--%%| <npVWRBHKNy|ukKeRjytSR>

# select a single array of data to use as a template
lwdown_data = ds_subset.isel(time=0).LWDOWN

# create a grid for the geocube
out_grid = make_geocube(
            vector_data=gdf,
            measurements=["cat"],
            like=ds_subset # ensure the data are on the same grid
)

# add the catchment variable to the original dataset
ds_subset = ds_subset.assign_coords(cat = (['y','x'], out_grid.cat.data))

# compute the unique catchment IDs which will be used to compute zonal statistics
catchment_ids = numpy.unique(ds_subset.cat.data[~numpy.isnan(ds_subset.cat.data)])

print(f'The dataset contains {len(catchment_ids)} catchments')

# |%%--%%| <ukKeRjytSR|2OsWn63isB>
r"""°°°
Plot the catchments of the subsetted area along with the grid cells corresponding to the subcatchments
°°°"""
# |%%--%%| <2OsWn63isB|OtMrAvpvET>

figure, ax = plt.subplots(figsize=(10,7))

# plot the gridded catchment mapping
ds_subset.cat.plot()

# preview map geometries
gdf.iloc[:].plot(ax=ax, linewidth=2, edgecolor='k', facecolor='None');

# |%%--%%| <OtMrAvpvET|3Y4il0AM1Z>
r"""°°°
Plot the catchments of the subsetted area along with the buffered cells
°°°"""
# |%%--%%| <3Y4il0AM1Z|D3yEzZkZUP>

figure, ax = plt.subplots(figsize=(10,7))

# plot cell_number to show the buffered area that has been selected
ds_subset.cell_number.plot()

# preview map geometries
gdf.iloc[:].plot(ax=ax, linewidth=2, edgecolor='k', facecolor='None');

# |%%--%%| <D3yEzZkZUP|CWLbgwgFUq>
r"""°°°
## Test using cell_number to compute weighted variables
°°°"""
# |%%--%%| <CWLbgwgFUq|Jk2TSuiu28>

# read the pre-computed weights 
weights = pandas.read_parquet('s3://lynker-spatial/v20.1/forcing_weights.parquet')

# |%%--%%| <Jk2TSuiu28|dYHL1MRe91>

# isolate a single catchment that we're interested in
cat_id = 2851705

# |%%--%%| <dYHL1MRe91|Pi59JfHycw>

# isolate a single catchment that we're interested in
cat_id = 2851705

# Select the cells that intersect with this boundary.
ds_catchment = ds.rio.clip(gdf.loc[gdf.divide_id == f'cat-{cat_id}'].geometry.values,
                 gdf.crs,
                 drop=True,
                 invert=False,
                 from_disk=True,
                 all_touched=True
                )

# |%%--%%| <Pi59JfHycw|VdyQBTJIpn>

ds_subset.where(ds_subset.cat==cat_id, drop=True).cell_number

# |%%--%%| <VdyQBTJIpn|iGRF9XGRkZ>

ds_catchment.cell_number

# |%%--%%| <iGRF9XGRkZ|ABxRN5Gs0c>

# todo: plot this catchment and the cells corresponding to the cell id.

figure, ax = plt.subplots(figsize=(10,7))


ds_catchment.cell_number.plot(levels=50)

# preview map geometries
current_geom = gdf.loc[gdf.divide_id == f'cat-{cat_id}']
gdf.loc[gdf.divide_id == f'cat-{cat_id}'].plot(ax=ax, linewidth=2, edgecolor='k', facecolor='None');

xmin = current_geom.geometry.bounds.minx.item()
xmax = current_geom.geometry.bounds.maxx.item()
ymin = current_geom.geometry.bounds.miny.item()
ymax = current_geom.geometry.bounds.maxy.item()

ax.set_ylim((ymin, ymax))
ax.set_xlim((xmin, xmax))

# |%%--%%| <ABxRN5Gs0c|L6VzoghCKU>

ds_catchment.cell_number

# |%%--%%| <L6VzoghCKU|XeXrQOTd3Z>

catchment_weights = weights.loc[weights.divide_id == f'cat-{cat_id}'].sort_values(by='cell')

ids = numpy.unique(ds_catchment.cell_number)
ids.sort()
for cell in ids:
    print(f'{cell} in weights -> {cell in catchment_weights.cell}')

# |%%--%%| <XeXrQOTd3Z|oTO4fYJdrm>

weights.loc[weights.divide_id == f'cat-{cat_id}'].sort_values(by='cell')

# |%%--%%| <oTO4fYJdrm|FPRkuPnoQH>



# |%%--%%| <FPRkuPnoQH|oFIj68DAjs>



# |%%--%%| <oFIj68DAjs|J1FHfjUyt0>

ymin.item()

# |%%--%%| <J1FHfjUyt0|b1tDZeaPfk>

# todo: plot this catchment and the cells corresponding to the cell id.

figure, ax = plt.subplots(figsize=(10,7))


# plot the gridded catchment mapping
#ds_catchment.cat.plot()

ds_catchment.cell_number.plot(levels=50)

# preview map geometries
current_geom = gdf.loc[gdf.divide_id == f'cat-{cat_id}']
gdf.loc[gdf.divide_id == f'cat-{cat_id}'].plot(ax=ax, linewidth=2, edgecolor='k', facecolor='None');

xmin = current_geom.geometry.bounds.minx.item()
xmax = current_geom.geometry.bounds.maxx.item()
ymin = current_geom.geometry.bounds.miny.item()
ymax = current_geom.geometry.bounds.maxy.item()

ax.set_ylim((ymin, ymax))
ax.set_xlim((xmin, xmax))

# |%%--%%| <b1tDZeaPfk|c1bxjGTPI7>

ds_catchment.cell_number.count()

# |%%--%%| <c1bxjGTPI7|xQXSDz3ogA>

weights.loc[weights.cell == 12337555]

# |%%--%%| <xQXSDz3ogA|DvpSi8sk9b>

idxmin = weights.cell.sub(12335252.2120625).abs().idxmin()
weights.loc[idxmin]

# |%%--%%| <DvpSi8sk9b|UDFJ89MJZS>

xs = ds_catchment.x.values

# |%%--%%| <UDFJ89MJZS|0scQ9MpdM8>

ys = ds_catchment.y.values

# |%%--%%| <0scQ9MpdM8|I1ghJZIReU>

y = ys - dy/2
x = xs - dx/2


# |%%--%%| <I1ghJZIReU|1isggqXsHJ>

(ymax - (y - dy/2)) / dy * ds_ncols + ((xs[0] - dx/2) - xmin) / dx

# |%%--%%| <1isggqXsHJ|Yqxg2PA35z>



# |%%--%%| <Yqxg2PA35z|RTuSw4aOo3>
r"""°°°
## Testing computing cell_values after subsetting region
°°°"""
# |%%--%%| <RTuSw4aOo3|3hA8KbQBEE>

# clip AORC to the extent of the hydrofabric geometries
ds_catchment = ds.rio.clip(gdf.geometry.values,
                 gdf.crs,
                 drop=True,
                 invert=False,
                 from_disk=True,
                 all_touched=True
                )

# |%%--%%| <3hA8KbQBEE|swnEHExdZe>

# isolate a single catchment that we're interested in
cat_id = 2851705
subcatchment = ds.where(ds.cat==cat_id, drop=True)

# |%%--%%| <swnEHExdZe|cR612S3kPj>

subcatchment

# |%%--%%| <cR612S3kPj|sPwq7PSuM9>

wgrid = json.load(open('forcing_grids.json', 'r'))

# |%%--%%| <sPwq7PSuM9|AgVKc5w8KK>

# only grab the x and y coordinates for a small region
xs = subcatchment.x.values
ys = subcatchment.y.values

# get the rest of the data from the large dataset
dx = ds.x.resolution
dy = ds.y.resolution
xmin = ds.x.min().item() - dx / 2 
xmax = ds.x.max().item() + dx / 2 
ymin = ds.y.min().item() - dy / 2 
ymax = ds.y.max().item() + dy / 2 
ds_nrows, ds_ncols = ds.LWDOWN.shape[1:]

# |%%--%%| <AgVKc5w8KK|cAwkfWgAct>

# unravel x, y coordinates into point pairs

rr, cc = numpy.meshgrid(xs, ys)
pair_pts = numpy.array([rr, cc]).T.reshape(-1, 2)

# |%%--%%| <cAwkfWgAct|qrg67APpte>

xss = pair_pts[:,0]
yss = pair_pts[:,1]

# |%%--%%| <qrg67APpte|UJGq9BSCEX>

#cell_number = ((ymax - (yss)) / dy * ds_ncols + ((xss) - xmin) / dx).round()
#cell_number = ((ymax - (yss - dy/2)) / dy * ds_ncols + ((xss- dy/2) - xmin) / dx).round()
cell_number = (ymax - (yss - dy/2)) / dy * ds_ncols + ((xss- dy/2) - xmin) / dx

# |%%--%%| <UJGq9BSCEX|IaicP3SWXe>

cell_ids = cell_number.reshape(ds_nrows, ds_ncols)
cell_ids.shape

# |%%--%%| <IaicP3SWXe|yJ7CKgmBKR>

cell_number

# |%%--%%| <yJ7CKgmBKR|xjAVOCwmA9>

cell_number

# |%%--%%| <xjAVOCwmA9|1OsmhtIl1V>

weights.loc[weights.divide_id == f'cat-{cat_id}']

# |%%--%%| <1OsmhtIl1V|6pQsIUdAdf>

weights.loc[weights.cell.isin(cell_number)]

# |%%--%%| <6pQsIUdAdf|zYQT3wOdQA>


