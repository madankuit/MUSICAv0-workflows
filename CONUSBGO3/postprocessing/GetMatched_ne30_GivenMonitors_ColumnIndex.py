'''
Extract the matched ne30np4 column index for a list of provided monitors.

For every monitor in the CONUSBGO3 monitor list, find the nearest MUSICA
ne30np4 model column and write the monitor -> column mapping to CSV.

All paths come from config/paths.py; nothing here is hard-coded.

MODIFICATION HISTORY:
    3 Sep 2025: VERSION 1.0
    - Initial version
    31 Aug 2026: VERSION 1.1
    - Paths moved to config/paths.py; SCRIP grid now read from the repo
    9 Oct 2026: VERSION 1.2
    - --experiment ne0CONUS matches the same monitors on the ne0CONUSne30x8 grid;
      prints the monitor-to-column-centre distance as a sanity check

    python GetMatched_ne30_GivenMonitors_ColumnIndex.py --experiment ne0CONUS
'''
#================================================================================================
# ### functions import ###

# basics
import os
import os.path
from os.path import exists
from os import path
import sys
import glob
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

#================================================================================================
# Configuration - every path is imported from config/paths.py
import pathlib
_ROOT = next(p for p in pathlib.Path(__file__).resolve().parents
             if (p / 'config' / 'paths.py').exists())
sys.path.insert(0, str(_ROOT))
import config  # noqa: F401  - also puts functions/ on sys.path
from config.paths import (
    BGO3_MONITOR_LIST,
    BGO3_EXPERIMENTS,
    case_hist_dir,
    ensure_dir,
)
import argparse
ap = argparse.ArgumentParser()
ap.add_argument('--experiment', default='ne30', choices=list(BGO3_EXPERIMENTS))
EXP = BGO3_EXPERIMENTS[ap.parse_args().experiment]

# repo-local shared functions (functions/ is on sys.path via `import config`)
from SE_analysis import get_site_index

# Specified input
MonitorInfo_filepath = BGO3_MONITOR_LIST
SCRIP_ne30 = str(EXP['scrip'])          # SE grid of this experiment (name kept from v1.0)
Savefile_path = EXP['colidx']
ensure_dir(Savefile_path.parent)

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
# Any h2 file of the experiment's first BASE case serves to read the lat/lon coordinates.
_ex_case = EXP['cases'][('BASE', EXP['years'][0])]
_ex_hist = case_hist_dir(_ex_case)
_ex_files = sorted(glob.glob(str(_ex_hist / '*.cam.h2.*.nc')))
if not _ex_files:
    raise FileNotFoundError(f'No h2 history files found under {_ex_hist}')
ex_h2_filepath = _ex_files[0]
ds_ne30 = xr.open_dataset(ex_h2_filepath).isel(lev=lev_idx,ilev=lev_idx,time=10) 
assert ds_ne30.sizes['ncol'] == EXP['ncol'], (ex_h2_filepath, ds_ne30.sizes['ncol'])

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
_ok = unique_monitor_locations['MUSICA0_colIndex'] != 'Find None'
_m = unique_monitor_locations[_ok]
_dlon = (((_m['Approx_MUSICA0_lon'].astype(float) - _m['lon'] + 180) % 360) - 180) * np.cos(np.deg2rad(_m['lat']))
_d_km = 111.2 * np.hypot(_m['Approx_MUSICA0_lat'].astype(float) - _m['lat'], _dlon)
print(f"{EXP['grid']}: matched {_ok.sum()}/{len(_ok)} monitors; monitor-to-column-centre distance "
      f"median {_d_km.median():.1f} km, p99 {_d_km.quantile(.99):.1f} km, max {_d_km.max():.1f} km")
unique_monitor_locations.to_csv(Savefile_path, index=False)

print("Saved to:",Savefile_path)

#================================================================================================
# 