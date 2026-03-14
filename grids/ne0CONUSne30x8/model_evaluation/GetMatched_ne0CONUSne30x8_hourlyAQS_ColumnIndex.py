'''
This script is used to extract matched column index in the ne0CONUSne30x8 horizontal grid of MUSICA model simulations matched with SLAMS/AQS monitors

MODIFICATION HISTORY:
    M. Tao, 18, DEC, 2023: VERSION 1.0
    - Initial version
'''
# ============================================================
# USER CONFIGURATION — set these paths before running
# ============================================================
AQS2018_diri = ''       # Path to AQS data directory for year 2018
AQS_diri = ''           # Path to base AQS data directory (contains parameters.csv)
MUSICA_index_diri = ''  # Path to directory for output column-index CSV files
SpinUp_diri = ''        # Path to CESM archive directory containing spinup h1 files
ex_h1_filename = ''     # Example h1 filename to read model lat/lon coordinates
# ============================================================

#================================================================================================
# Specified variables
hourly_varlist = ['CO','SO2','NO2','O3','LC25']
# 'LC25' for PM2.5 - Local Conditions (88101)

# for most hourly measurement
valcolname = 'Sample Measurement'

#================================================================================================
# ### functions import ###

# basics
import os
import os.path
from os.path import exists
from os import path
import sys  
from io import BytesIO
import requests
from zipfile import ZipFile
import json

# numerical
import math
import random
import numpy as np
import numpy.ma as ma
import pandas as pd
import netCDF4 as nc4
from netCDF4 import Dataset
import xarray as xr

# my functions
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../functions/'))
from func_ModelEval_statistical_tests import *
from SE_analysis import get_site_index

# Grid
# Read SCRIP file that has grid information needed to plot values on a map
SCRIP_CONUS = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../grid_files/ne0CONUS_ne30x8_np4_SCRIP.nc')
SCRIP_ne30 = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../grid_files/ne30np4_091226_pentagons.nc')

#================================================================================================
### """Functions for AQS"""
from zipfile import ZipFile

def deflatten(names):
    """This function create a sort of structure from the file list of a zip file"""
    names.sort(key=lambda name:len(name.split('/')))
    deflattened = []
    while len(names) > 0:
        name = names[0]
        if name[-1] == '/':
            subnames = [subname[len(name):] for subname in names if subname.startswith(name) and subname != name]
            for subname in subnames:
                names.remove(name+subname)
            deflattened.append((name, deflatten(subnames)))
        else:
            deflattened.append(name)
        names.remove(name)
    return deflattened

# get parameter code
parameters_df = pd.read_csv(AQS_diri+"parameters.csv")
ParaAbrrCode_dic = dict(zip(parameters_df["Parameter Abbreviation"].values, parameters_df["Parameter Code"].values))
ParaAbrrName_dic = dict(zip(parameters_df["Parameter Abbreviation"].values, parameters_df["Parameter"].values))
ParaAbrrUnit_dic = dict(zip(parameters_df["Parameter Abbreviation"].values, parameters_df["Standard Units"].values))

# "Formaldehyde" in parameters_df["Parameter"].unique()
# HCHO_code = ParaNameCode_dic["Formaldehyde"]
# ISOP_code = ParaNameCode_dic["Isoprene"]    
#================================================================================================
# Function to get the MUSICA column index for hourly dataframe
def write_hourlyAQSMonitorID_to_MUSICAcolidx_csvfile(varname,AQSvar_df,CVSOUT_diri):
    """
    The function calculates the approximate column index of the MUSICA grid for each unique MonitorID location and saves the results to a CSV file.
    Inputs:
        varname (str): The name of the species or variable for which column indices are being calculated.
        AQSvar_df (pandas.DataFrame): A DataFrame containing data for the specified species, including 'MonitorID', 'Latitude', and 'Longitude'.
        CVSOUT_diri: a path to a directory to save

    CSVs are saved in the 'colidxCSV' directory with filenames like MUSICA_index_diri+'MUSICA0_hourlyAQS_'+{varname}+'_Monitor_colidx.csv'
    Each CSV contains columns:
    - 'MonitorID': Unique MonitorID values
    - 'Latitude': Latitude values for MonitorID locations
    - 'Longitude': Longitude values for MonitorID locations
    - 'MUSICA0_colIndex': Approximate column index on the MUSICA grid
    - 'Approx_MUSICA0_lat': Approximate latitude on the MUSICA grid
    - 'Approx_MUSICA0_lon': Approximate longitude on the MUSICA grid
    """
    
    # check if the file already existed
    Savefile_path = CVSOUT_diri+'MUSICA0_hourlyAQS_'+varname+'_Monitor_colidx.csv'
    if os.path.exists(Savefile_path):
        print(f"Skip processing {varname}")
    else:
        print(f"Start....")
        # Read an example MUSICA ne0CONUSne30x8 ds for model lat and lon
        ds_CONUS = xr.open_dataset(SpinUp_diri+ex_h1_filename)

        # Keep only the unique MonitorID values along with their Latitude and Longitude
        unique_AQSmonitor_locations = AQSvar_df[['MonitorID', 'Latitude', 'Longitude']].drop_duplicates().reset_index(drop=True)

        # loop through the monitors and find the corresponding column index
        ls_Index_MonitorIDi = []
        ls_Model_lati = []
        ls_Model_loni = []

        for MonitorIDi in unique_AQSmonitor_locations.MonitorID.values:
            # Get the latitude and longitude values for 'MonitorIDi'
            MonitorIDi_df = unique_AQSmonitor_locations[unique_AQSmonitor_locations['MonitorID'] == MonitorIDi]

            lati = MonitorIDi_df['Latitude'].iloc[0]
            loni = MonitorIDi_df['Longitude'].iloc[0]

            # find the model index
            Index_MonitorIDi = get_site_index( site_lat=lati, site_lon=360+loni, scrip_file=SCRIP_CONUS )
            if Index_MonitorIDi==None:
                # add to the list; use np.nan as sentinel since Model_lati/loni are not yet defined
                ls_Index_MonitorIDi.append('Find None')
                ls_Model_lati.append(np.nan)
                ls_Model_loni.append(np.nan)

            else:
                Model_lati = ds_CONUS.lat.values[Index_MonitorIDi]
                Model_loni = ds_CONUS.lon.values[Index_MonitorIDi]
                # add to the list
                ls_Index_MonitorIDi.append(Index_MonitorIDi)
                ls_Model_lati.append(Model_lati)
                ls_Model_loni.append(Model_loni)

        # append to the df
        unique_AQSmonitor_locations['MUSICA0_colIndex'] = ls_Index_MonitorIDi
        unique_AQSmonitor_locations['Approx_MUSICA0_lat'] = ls_Model_lati
        unique_AQSmonitor_locations['Approx_MUSICA0_lon'] = ls_Model_loni

        # Save the DataFrame to a CSV file
        unique_AQSmonitor_locations.to_csv(Savefile_path, index=False)

        print("Saved to:",Savefile_path)
    
#================================================================================================
### Read in hourly concentrations for given variables and find the column index in MUSICA ne0CONUSne30x8
for varname in hourly_varlist:
    print(varname)
    #-------------------------------------------------------
    # Get AQS
    AQS_filename = 'hourly_'+str(ParaAbrrCode_dic[varname])+'_2018.zip'
    AQS_zf = ZipFile(AQS2018_diri+AQS_filename)
    # see what's in the zip file
    Inside_filenames = deflatten(AQS_zf.namelist())
    # read in the csv file in the zip file ,
    fields = ['State Code','County Code','Site Num','Parameter Code','POC','Latitude','Longitude','Date Local',valcolname,'Units of Measure']
    AQS_df = pd.read_csv(AQS_zf.open(Inside_filenames[0]),skipinitialspace=True, index_col=False, usecols=fields)
    # # Add a new column combining State Code-County Code-Site Num
    AQS_df = AQS_df.astype({col: str for col in ["State Code", "County Code", "Site Num", "Parameter Code"]})
    AQS_df["MonitorID"] = AQS_df[["State Code", "County Code", "Site Num","Parameter Code"]].agg('-'.join, axis=1) 
    # Simplify the df
    AQS_df = AQS_df[['MonitorID','Latitude','Longitude','Date Local',valcolname,'Units of Measure']]
    # Keep only the unique MonitorID values along with their Latitude and Longitude
    unique_AQSmonitor_locations = AQS_df[['MonitorID', 'Latitude', 'Longitude']].drop_duplicates().reset_index(drop=True)
    #-------------------------------------------------------
    # Find the matched index
    write_hourlyAQSMonitorID_to_MUSICAcolidx_csvfile(varname,unique_AQSmonitor_locations,MUSICA_index_diri)

#================================================================================================
# 