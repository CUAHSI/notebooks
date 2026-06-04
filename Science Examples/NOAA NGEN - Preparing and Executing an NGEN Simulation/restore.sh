#!/bin/bash

while true; do

echo "
This script will revert changes made to the hydrofabric data 
located in ~/.ngiab/hydrofabric/v2.2/.
"
    
    read -p "Do you want to proceed? (y/n): " yn
    case $yn in
        [Yy]* ) break;;
        [Nn]* ) exit;;
        * ) continue ;;
    esac
done

# Remove the files we copied into this directory
echo "Removing VPU 03W hydrofabric files"
if [ -f "/home/jovyan/.ngiab/hydrofabric/v2.2/conus_nextgen.gpkg" ]; then
  rm "/home/jovyan/.ngiab/hydrofabric/v2.2/conus_nextgen.gpkg" 
fi
if [ -f "/home/jovyan/.ngiab/hydrofabric/v2.2/conus_igraph_network.gpickle" ]; then
  rm "/home/jovyan/.ngiab/hydrofabric/v2.2/conus_igraph_network.gpickle" 
fi

# Rename the original hydrofabric files which 
# should have a .bak suffix
echo "Renaming existing hydrofabric files"
if [ -f "/home/jovyan/.ngiab/hydrofabric/v2.2/conus_nextgen.gpkg.bak" ]; then
  mv "/home/jovyan/.ngiab/hydrofabric/v2.2/conus_nextgen.gpkg.bak" "/home/jovyan/.ngiab/hydrofabric/v2.2/conus_nextgen.gpkg"
fi

if [ -f "/home/jovyan/.ngiab/hydrofabric/v2.2/conus_igraph_network.gpickle.bak" ]; then
  mv "/home/jovyan/.ngiab/hydrofabric/v2.2/conus_igraph_network.gpickle.bak" "/home/jovyan/.ngiab/hydrofabric/v2.2/conus_igraph_network.gpickle"
fi
