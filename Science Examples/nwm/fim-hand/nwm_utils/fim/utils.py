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
from typing import Dict, List
from pathlib import Path
from scipy import interpolate
from matplotlib import pyplot as plt


class Reach:
    def __init__(self, reach_id, flows=[], times=[], hydroids=[]):
        self._id = reach_id
        self._flows = flows
        self._times = times
        self._hydroids = hydroids
        self._rating_curve_dict = {}

    # accessor
    def id(self):
        return self._id
    def flows(self, flows=None):
        if flows is None:
            return self._flows
        self._flows = flows
    def times(self, times=None):
        if times is None:
            return self._times
        self._times = times
    def hydroids(self, hydroids=None):
        if hydroids is None:
            return self._hydroids
        self._hydroids = hydroids
    def get_rating_curve(self, reach_id=None):
        if reach_id is None:
            return self._rating_curve_dict
        return self._rating_curve_dict[reach_id]

    # publics
    def assign_flows(self, flows, times):
        if len(flows) != len(times):
            raise Exception('len()flows) != len(times)')
        self.times(times)
        self.flows(flows)
    def set_rating_curve_from_hydrotable(self, hydrotable):
        # get unique combos of HydroID and NextDownID
        # set these in the class variable
        df = hydrotable.loc[hydrotable.feature_id == self._id]
        self.hydroids(np.unique(df.HydroID))

        # collect rating curve data for each hydroid
        for hydro_id in self.hydroids():
            print(f'Setting RC for {hydro_id}')
            # look up rating curve for this hydroid
            dat = df.loc[df.HydroID == hydro_id, ["stage", "discharge_cms"]]
            rc = RatingCurve(hydro_id,
                             stage = dat.stage,
                             discharge = dat.discharge_cms)
            self._rating_curve_dict[hydro_id] = rc

class RatingCurve:
    def __init__(self, reach_id, stage, discharge):
        self._id = reach_id
        self._stage = stage
        self._discharge = discharge
        self._interp_func = self.__prepare_interp_func()

    # accessors
    def stage(self, stage=None):
        if stage is None:
            return self._stage
        self._stage = stage
    def discharge(self, discharge=None):
        if discharge is None:
            return self._discharge
        self._discharge = discharge

    # privates
    def __prepare_interp_func(self):
        return interpolate.interp1d(self.discharge(),
                                 self.stage(), kind="linear")

    # publics
    def plot(self):
        fig, ax = plt.subplots()
        ax.plot(self.discharge(), self.stage());
        ax.set_xlabel('Discharge [cms]')
        ax.set_ylabel('Stage [m]')
        ax.set_title(f'Rating Curve for {self._id}')
        ax.grid(True)
        fig.tight_layout()

    def get_stage_for_flow(self, q):
        return self._interp_func(q)


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


#def get_stage_for_nhd_reach(
#    nhd_feature_id: int, flow_cms: float, hydrotable: Path = Path("./hydroTable_0.csv")
#) -> Dict[int, float]:
#    """
#    Retrieves stage for a given NHD reach and input streamflow. The stage for
#    the most downstream HydroID is computed using FIM rating curves.
#
#    Parameters
#    ==========
#    nhd_feature_id: int
#        NHD feature identifier.
#    flow_cms: float
#        Streamflow for the reach in cubic meters per second
#    hydrotable: pathlib.Path
#        Path to the FIM hydrotable.csv file containing rating curve data.
#
#    Returns
#    =======
#    Dict [int, float]
#        Dictionary containing a single NWM hydro identifier and
#        its corresponding stage.
#
#    """
#
#    # load hydrotable_0
#    # we don't need all of the columns in this csv
#    hydro_df = pd.read_csv(
#        hydrotable,
#        usecols=["HydroID", "NextDownID", "feature_id", "stage", "discharge_cms"],
#    )
#
#    # select features that match nhd_feature_id
#    d = hydro_df.loc[hydro_df.feature_id == nhd_feature_id]
#
#    # get unique combos of HydroID and NextDownID
#    from_nodes, to_nodes = np.unique(d[["HydroID", "NextDownID"]], axis=0).T
#
#    # find the most downstream reach
#    # computed as the 'From_Node' corresponding to the 'To_Node' that
#    # that does not exist in the list of 'From_Node'
#    downstream_hydroid = None
#    for i, to_node in enumerate(to_nodes):
#        if to_node not in from_nodes:
#            # this is a node beyond our NHD reach
#            # save the corresponding "from_node" id
#            downstream_hydroid = from_nodes[i]
#
#    # TODO: insert error handling here
#    if downstream_hydroid == None:
#        print("Something went wrong :( ")
#        return {-9999: -9999.0}
#
#    interpolated_stage = compute_stage(d, downstream_hydroid, flow_cms)
#
#    return interpolated_stage


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

    # # load hydrotable_0
    # # we don't need all of the columns in this csv
    # hydro_df = pd.read_csv(
    #     hydrotable,
    #     usecols=["HydroID", "NextDownID", "feature_id", "stage", "discharge_cms"],
    # )

    # # select features that match nhd_feature_id
    # d = hydro_df.loc[hydro_df.feature_id == nhd_feature_id]

    
    d = get_hydroids_for_reach(nhd_feature_id, hydrotable)

    # get unique combos of HydroID and NextDownID
    hydro_ids = np.unique(d.HydroID)

    interpolated_stages = {}
    for hydro_id in hydro_ids:
        interpolated_stages.update(compute_stage(d, hydro_id, flow_cms))

    # return interpolated stage
    return interpolated_stages

def get_hydroids_for_reach(nhd_feature_id: int, hydrotable: Path = Path("./hydroTable_0.csv")) -> pd.DataFrame:

    # load hydrotable_0
    # we don't need all of the columns in this csv
    hydro_df = pd.read_csv(
        hydrotable,
        usecols=["HydroID", "NextDownID", "feature_id", "stage", "discharge_cms"],
    )

    # select features that match nhd_feature_id
    d = hydro_df.loc[hydro_df.feature_id == nhd_feature_id]

    return d
    # # get unique combos of HydroID and NextDownID
    # hydro_ids = np.unique(d.HydroID)

    

    return hydro_ids
