#!/usr/bin/env python3

"""
Purpose: Utility functions for working with NOAA OWP FIM data
Authors: Tony Castronova <acastronova@cuahsi.org>
         Irene Garousi-Nejad <igarousi@cuahsi.org>
Last Modified: March 26, 2024
"""


import xarray
import rioxarray
import numpy as np
import pandas as pd
from typing import Dict
from pathlib import Path
from scipy import interpolate

def get_hand_object(path: Path) -> xarray.Dataset:
    """
    Converts NOAA OWP HAND tif into a rioxarray Dataset.

    Parameters
    ==========
    path: pathlib.Path
        Path-like object to the NOAA OWP HAND tif.

    Returns
    =======
    xarray.Dataset
        Dataset containing HAND as a variable.

    """
    return rioxarray.open_rasterio(path, masked=True) \
            .squeeze() \
            .drop_vars('band'). \
            to_dataset(name='hand')

def get_nwm_id(geo_lat, geo_lon, ds):
    """
    This function determines the NWM ID linked to a nearby reach based on a specified geographic location.
    input:
        geo_lat: latitude of the geo-location
        geo_lon: longitude of the geo-location
        ds: the NWM dataset that includes streamflow data
    output:
        NWM reach/feature ID
    """
    # Calculate the Euclidean distance to find the closest point
    distances = np.sqrt((ds.latitude - geo_lat) ** 2 + (ds.longitude - geo_lon) ** 2)

    # Find the index of the minimum distance
    min_distance_index = np.unravel_index(distances.argmin(), distances.shape)

    # Extract the subset based on the indices
    subset = ds.isel(feature_id=min_distance_index[0].item()).compute()

    return subset.feature_id.values


def interpolate_y(
    df: pd.DataFrame, x_column: str, y_column: str, x_value: float
) -> float:
    """
    Performs 1D interpolation on two columns of a Dataframe.

    Parameters
    ==========
    df: pandas.DataFrame
        DataFrame containing data that will be used in the interpolation.
    x_column: str
        Name of the column that represents the X-axis data.
    y_column: str
        Name of the column that represents the Y-axis data.
    x_value: float
        Numeric X-axis value for which to interpolate Y-axis data. Returns
        -9999 if the interpolation fails to resolve.

    Returns
    =======
    y_value: float
        Numeric Y-axis value corresponding to the input X-axis value.

    """
    # Sort the DataFrame by the 'x' column to ensure interpolation works correctly
    df_sorted = df.sort_values(by=x_column)

    # Check if the x_value is within the range of the DataFrame
    if x_value < df_sorted[x_column].min() or x_value > df_sorted[x_column].max():
        return -9999  # x_value is out of range, cannot interpolate

    # Perform linear interpolation
    f = interpolate.interp1d(df_sorted[x_column], df_sorted[y_column], kind="linear")

    # Return the interpolated y value for the given x value
    return f(x_value)


def compute_stage(df: pd.DataFrame, hydro_id: int, flow_cms: float) -> Dict[int, float]:
    """
    Computes river stage from a rating curve given streamflow in cms.

    Parameters
    ==========
    df: pandas.DataFrame
        DataFrame containing the stage and discharge values of the rating curve.
        This must contain the following columns: HydroID, stage, discharge_cms.
    hydro_id: int
        Identifier for the reach for which to compute stage.
    flow_cms: float
        Streamflow to convert into river stage.

    Returns
    =======
    Dict [int, float]
        A dictionary containing computed stage and its associated hydroid

    """

    # look up rating curve for this hydroid
    rating_curve = df.loc[df.HydroID == hydro_id, ["stage", "discharge_cms"]]

    # interpolate using the provided flow rate
    interpolated_stage = interpolate_y(rating_curve, "discharge_cms", "stage", flow_cms)

    return {hydro_id: float(interpolated_stage)}


def get_stage_for_nhd_reach(
    nhd_feature_id: int, flow_cms: float, hydrotable: Path = Path("./hydroTable_0.csv")
) -> Dict[int, float]:
    """
    Retrieves stage for a given NHD reach and input streamflow. The stage for
    the most downstream HydroID is computed using FIM rating curves.

    Parameters
    ==========
    nhd_feature_id: int
        NHD feature identifier.
    flow_cms: float
        Streamflow for the reach in cubic meters per second
    hydrotable: pathlib.Path
        Path to the FIM hydrotable.csv file containing rating curve data.

    Returns
    =======
    Dict [int, float]
        Dictionary containing a single NWM hydro identifier and
        its corresponding stage.

    """

    # load hydrotable_0
    # we don't need all of the columns in this csv
    hydro_df = pd.read_csv(
        hydrotable,
        usecols=["HydroID", "NextDownID", "feature_id", "stage", "discharge_cms"],
    )

    # select features that match nhd_feature_id
    d = hydro_df.loc[hydro_df.feature_id == nhd_feature_id]

    # get unique combos of HydroID and NextDownID
    from_nodes, to_nodes = np.unique(d[["HydroID", "NextDownID"]], axis=0).T

    # find the most downstream reach
    # computed as the 'From_Node' corresponding to the 'To_Node' that
    # that does not exist in the list of 'From_Node'
    downstream_hydroid = None
    for i, to_node in enumerate(to_nodes):
        if to_node not in from_nodes:
            # this is a node beyond our NHD reach
            # save the corresponding "from_node" id
            downstream_hydroid = from_nodes[i]

    # TODO: insert error handling here
    if downstream_hydroid == None:
        print("Something went wrong :( ")
        return {-9999: -9999.0}

    interpolated_stage = compute_stage(d, downstream_hydroid, flow_cms)

    return interpolated_stage


def get_stage_for_all_hydroids_in_reach(
    nhd_feature_id: int, flow_cms: float, hydrotable: Path = Path("./hydroTable_0.csv")
) -> Dict[int, float]:
    """
    Retrieves stage for all NWM HydroIDs given an NHD reach and input streamflow.
    The stage is computed using FIM rating curves.

    Parameters
    ==========
    nhd_feature_id: int
        NHD feature identifier.
    flow_cms: float
        Streamflow for the reach in cubic meters per second
    hydrotable: pathlib.Path
        Path to the FIM hydrotable.csv file containing rating curve data.

    Returns
    =======
    Dict [int, float]
        Dictionary containing one or more NWM hydro identifiers and
        their corresponding stages.

    """

    # load hydrotable_0
    # we don't need all of the columns in this csv
    hydro_df = pd.read_csv(
        hydrotable,
        usecols=["HydroID", "NextDownID", "feature_id", "stage", "discharge_cms"],
    )

    # select features that match nhd_feature_id
    d = hydro_df.loc[hydro_df.feature_id == nhd_feature_id]

    # get unique combos of HydroID and NextDownID
    hydro_ids = np.unique(d.HydroID)

    interpolated_stages = {}
    for hydro_id in hydro_ids:
        interpolated_stages.update(compute_stage(d, hydro_id, flow_cms))

    # return interpolated stage
    return interpolated_stages
