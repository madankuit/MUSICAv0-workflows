'''
This script is used to extract matched column index in the ne30 horizontal grid of MUSICA model simulations matched with a list of provided monitors

MODIFICATION HISTORY:
    Madankui Tao, 3, Sep, 2025: VERSION 1.0
    - Initial version
'''
#================================================================================================
# Specified input
Output_diri = '/net/fs09/d0/taoma528/ProcessedData/DanJaffeMUSICAPostprocessing/'

MonitorInfo_filepath = f'{Output_diri}MonitorInfo/Lee_Jaffe_GAM_stats.csv'

# Grid
SCRIP_ne30 = '/home/taoma528/Scripts/CESM_analysis/functions/ne30np4_091226_pentagons.nc'

Savefile_path = f'{Output_diri}MonitorInfo/MatchedMonitors_ne30_ColIdx.csv'

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
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../functions/'))  # repo-local functions/
from func_ModelEval_statistical_tests import *
from SE_analysis import get_site_index

#================================================================================================
### """Monitor Info"""
import pandas as pd
from pathlib import Path

def load_monitor_info(MonitorInfo_filepath) -> pd.DataFrame:
    """Load monitor metadata with columns: lon, lat, site_name, AQS_code."""
    path = Path(MonitorInfo_filepath)
    usecols = ["lon", "lat", "site_name", "AQS_code"]
    dtypes  = {"site_name": "string", "AQS_code": "string"}  # keep leading zeros

    if path.suffix.lower() in {".xlsx", ".xls"}:
        df = pd.read_excel(path, usecols=usecols, dtype=dtypes)
    else:
        # sep=None lets pandas infer comma/space/tab; engine='python' required for that
        df = pd.read_csv(path, usecols=usecols, dtype=dtypes, sep=None, engine="python")

    # Basic cleaning
    df.columns = [c.strip() for c in df.columns]
    df["site_name"] = df["site_name"].str.strip()
    df["AQS_code"]  = df["AQS_code"].str.strip()

    # Ensure numeric coords
    df["lat"] = pd.to_numeric(df["lat"], errors="coerce")
    df["lon"] = pd.to_numeric(df["lon"], errors="coerce")

    # Convert any 0–360 longitudes to −180–180
    needs_wrap = df["lon"] > 180
    if needs_wrap.any():
        df.loc[needs_wrap, "lon"] = df.loc[needs_wrap, "lon"] - 360

    # Drop rows missing coordinates
    df = df.dropna(subset=["lat", "lon"]).reset_index(drop=True)

    # Add a stable site index (0..N-1) you can carry into NetCDF as the 'site' dimension
    df.insert(0, "site", range(len(df)))

    return df

# Get the monitor info
MonitorInfo_df = load_monitor_info(MonitorInfo_filepath)
# Keep only the unique MonitorID values along with their Latitude and Longitude
unique_monitor_locations = MonitorInfo_df[['AQS_code', 'lat', 'lon']].drop_duplicates().reset_index(drop=True)

#================================================================================================
## Get an example MUSICA output file
lev_idx = -1 # for surface
ex_h2_filepath = '/net/fs09/d0/taoma528/CESM22/archive/f.e22.FCnudged.ne30_ne30_mg17.BGO3.BASEY20220401TY20230401/atm/hist/f.e22.FCnudged.ne30_ne30_mg17.BGO3.BASEY20220401TY20230401.cam.h2.2022-05-04-03600.nc'
ds_ne30 = xr.open_dataset(ex_h2_filepath).isel(lev=lev_idx,ilev=lev_idx,time=10) 

#================================================================================================
# Get the MUSICA column index for hourly dataframe

# loop through the monitors and find the corresponding column index
ls_Index_MonitorIDi = []
ls_Model_lati = []
ls_Model_loni = []

# # try one monitor or several
# for MonitorIDi in unique_monitor_locations.AQS_code.values[:5]:
# if all
for MonitorIDi in unique_monitor_locations.AQS_code.values:
    # Get the latitude and longitude values for 'MonitorIDi'
    MonitorIDi_df = unique_monitor_locations[unique_monitor_locations['AQS_code'] == MonitorIDi]

    lati = MonitorIDi_df['lat'].iloc[0]
    loni = MonitorIDi_df['lon'].iloc[0]

    # find the model index
    Index_MonitorIDi = get_site_index( site_lat=lati, site_lon=360+loni, scrip_file=SCRIP_ne30 )
    if Index_MonitorIDi==None:
        # add to the list
        ls_Index_MonitorIDi.append('Find None')
        ls_Model_lati.append(Model_lati)
        ls_Model_loni.append(Model_loni)

    else:
        Model_lati = ds_ne30.lat.values[Index_MonitorIDi]
        Model_loni = ds_ne30.lon.values[Index_MonitorIDi]
        # add to the list
        ls_Index_MonitorIDi.append(Index_MonitorIDi)
        ls_Model_lati.append(Model_lati)
        ls_Model_loni.append(Model_loni)
        
# append to the df
unique_monitor_locations['MUSICA0_colIndex'] = ls_Index_MonitorIDi
unique_monitor_locations['Approx_MUSICA0_lat'] = ls_Model_lati
unique_monitor_locations['Approx_MUSICA0_lon'] = ls_Model_loni

# Save the DataFrame to a CSV file
unique_monitor_locations.to_csv(Savefile_path, index=False)

print("Saved to:",Savefile_path)

#================================================================================================
# 