#!/bin/bash

while true; do

echo "
This script will place the hydrofabric data located in
~/.ngiab/hydrofabric/v2.2/, with a subset that is provided
in this HydroShare resource. This operation can be reverted
using the ./restore.sh utility.
"
    
    read -p "Do you want to proceed? (y/n): " yn
    case $yn in
        [Yy]* ) break;;
        [Nn]* ) exit;;
        * ) continue ;;
    esac
done

# make ngiab directory if it doesn't exist
mkdir -p "/home/jovyan/.ngiab/hydrofabric/v2.2/"

# Rename hydrofabric files if they exist
echo "Renaming existing hydrofabric files"
if [ -f "/home/jovyan/.ngiab/hydrofabric/v2.2/conus_nextgen.gpkg" ]; then
  mv "/home/jovyan/.ngiab/hydrofabric/v2.2/conus_nextgen.gpkg" "/home/jovyan/.ngiab/hydrofabric/v2.2/conus_nextgen.gpkg.bak"
fi

if [ -f "/home/jovyan/.ngiab/hydrofabric/v2.2/conus_igraph_network.gpickle" ]; then
  mv "/home/jovyan/.ngiab/hydrofabric/v2.2/conus_igraph_network.gpickle" "/home/jovyan/.ngiab/hydrofabric/v2.2/conus_igraph_network.gpickle.bak"
fi

# Copy the Hydrofabric subset from the resource into the ~/.ngiab folder. 
# The purpose of this is to eliminate the need to download 5+ GB of data
echo "Copying resource-specific hydrofabric files to ~/.ngiab/hydrofabric/v2.2/"
cp hydrofabric_vpu_03w/vpu-03W_subset.gpkg /home/jovyan/.ngiab/hydrofabric/v2.2/conus_nextgen.gpkg
cp hydrofabric_vpu_03w/conus_igraph_network.gpickle  /home/jovyan/.ngiab/hydrofabric/v2.2/conus_igraph_network.gpickle

echo "

Run the following command to see the files that have been changed:

ls -lah /home/jovyan/.ngiab/hydrofabric/v2.2

"