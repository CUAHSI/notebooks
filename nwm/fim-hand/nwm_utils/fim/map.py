#!/usr/bin/env python3


import time
import geopandas
import ipyleaflet
# import geopandas as gpd
from pathlib import Path
from sidecar import Sidecar
from ipywidgets import Layout

import shapely
import pynhd as nhd

import warnings
warnings.filterwarnings("ignore")

class SideCarMap():
    def __init__(self, basemap=ipyleaflet.basemaps.OpenStreetMap.Mapnik, gdf=None):
        self.selected_id = None
        self.map = None
        self.basemap = basemap
        self.gdf = gdf

    def display_map(self):
        defaultLayout=Layout(width='960px', height='940px')
        center = (45.9163, -94.8593)

        print('Creating Base Map...', end='')
        st = time.time()
        m = ipyleaflet.Map(
        basemap=ipyleaflet.basemap_to_tiles(ipyleaflet.basemaps.OpenStreetMap.Mapnik, layout=defaultLayout),
            center=center,
            zoom=9,
            scroll_wheel_zoom=True,
            tap=False
            )
        print(f'{time.time() - st:0.2f} sec')
        
        
        # add USGS Gages
        print('Adding USGS Gages...', end='')
        st = time.time()
        m.add_layer(
            ipyleaflet.WMSLayer(
                url='http://arcgis.cuahsi.org/arcgis/services/NHD/usgs_gages/MapServer/WmsServer',
                layers='0',
                transparent=True,
                format='image/png',
                min_zoom=8,
                max_zoom=18,
                )
        )
        print(f'{time.time() - st:0.2f} sec')
        
        # add NHD+ Reaches
        print('Adding NHD+ Reaches...', end='')
        st = time.time()
        m.add_layer(
            ipyleaflet.WMSLayer(
                url='https://hydro.nationalmap.gov/arcgis/services/nhd/MapServer/WMSServer',
                layers='6',
                transparent=True,
                format='image/png',
                min_zoom=8,
                max_zoom=18,
                )
        )
        print(f'{time.time() - st:0.2f} sec')

        # add features from geopandas if they are provided
        if self.gdf is not None:
            print('Adding Geopandas Features...', end='')
            st = time.time()
            approx_center = self.gdf.dissolve().centroid
            center=(approx_center.iloc[0].y, approx_center.iloc[0].x)
            
            geo_data = ipyleaflet.GeoData(geo_dataframe = self.gdf,
                   style={'color': 'blue', 'opacity':0.5, 'weight':1.9,}
                  )
            m.add(geo_data)

            # update the map center point
            m.center = center

            print(f'{time.time() - st:0.2f} sec')

        # bind the map handler function
        m.on_interaction(self.handle_map_interaction)

        print('Rendering Map...', end='')
        st = time.time()
        sc = Sidecar(title='NHD+ River Reaches')
        with sc:
            display(m)
        print(f'{time.time() - st:0.2f} sec')
        
        # save the map object for later
        self.map = m
        
    def handle_map_interaction(self, **kwargs):
    
        if kwargs.get('type') == 'click':
            print(kwargs)
            lat, lon = kwargs['coordinates']
            print(f'{lat}, {lon}')
            
            # query the reach nearest this point
            point = shapely.Point(lon, lat)

            # buffer the selected point by a small degree. This
            # is a hack for now and Buffer operations should only
            # be applied in a projected coordinate system in the future.
            print('buffering')
            pt_buf = point.buffer(0.001) 
            # wlayer = ipyleaflet.WKTLayer(
            #         wkt_string=pt_buf.wkt,
            #         )
            # self.map.add(wlayer);

            try:
                # remove the previously selected layers
                while len(self.map.layers) > 4:
                    self.map.remove(self.map.layers[-1]);
                
                # query the FIM reach that intersects with the point
                print('intersecting...')
                mask = self.gdf.intersects(pt_buf)
                print(f'found {len(self.gdf.loc[mask])} reaches')
                print('saving selected...')
                self.selected(value=self.gdf.loc[mask].iloc[0])

                # highlight this layer on the map
                wlayer = ipyleaflet.WKTLayer(
                    wkt_string=self.selected().geometry.wkt,
                    style={'color': 'green', 'opacity':1, 'weight':2.,})
                self.map.add(wlayer)
                
            except Exception: 
                print('Could not find reach for selected area')

    # getter/setter for the selected reach
    def selected(self, value=None):
        if value is None:
            if self.selected_id is None:
                print('No reach is selected.\nUse the map interface to select a reach of interest')
            else:
                return self.selected_id
        else:
            self.selected_id = value