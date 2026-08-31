#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extract_givenmonitorO3_hourly_dailyMDA8_toNetCDF.py

Description
-----------
Extract hourly surface O3 from MUSICA outputs at specified monitor locations,
compute daily maximum 8-hour average O3 (MDA8) following EPA convention, and
save both hourly (UTC) and daily (local-time) results into CF-compliant
NetCDF files. Reference to Extract_MDA8O3_givenMonitors.ipynb

Inputs
------
- casename: str
    Simulation case key for MUSICA outputs.
- MonitorIdx_df: pandas.DataFrame
    Site metadata table with columns:
      * lon, lat: coordinates of monitor
      * AQS_code: site identifier
      * MUSICA0_colIndex: column index in model output
- startMMDD, endMMDD: str
    Start and end dates in 'MM-DD' format (inclusive).
- Output_diri: str
    Directory where output NetCDF files are written.
- casename_label_dic: dict
    Maps casename to string label (last 4 chars = year).
- mergedh2_surfO3filepath_dic: dict
    Maps casename to merged hourly surface O3 file path.
- compute_mda8_UTCoffset: function
    Computes daily MDA8 O3 from UTC hourly series given UTC offset (hours).
- get_summer_offset: function
    Returns integer UTC offset hours for given lat/lon (fixed summertime DST).

Outputs
-------
Two NetCDF files in Output_diri:
- LocalTimeMDA8O3.<case_label>.GivenMonitors.<startDate>T<endDate>.nc
    Daily MDA8 O3 at all monitor sites (units: ppb).
- UTChourlyO3.<case_label>.GivenMonitors.<startDate>T<endDate>.nc
    Hourly surface O3 at all monitor sites (units: model original).

Notes
-----
- MDA8 calculation: rolling 8-h mean with ≥6 valid hours; daily max from
  windows ending 07–23 local; day valid if ≥13 valid windows.
- Output NetCDFs include per-site coordinates (lat, lon, AQS_code).
- Script prints the saved file paths; no return object.

MODIFICATION HISTORY:
    4, Sep, 2025: VERSION 1.0
    - Initial version
"""
# ============================================================
# CONFIGURATION - every path comes from config/paths.py, the single
# source of truth. Override cluster locations with the MUSICA_ENV_*
# environment variables documented there. Do not hard-code paths here.
# ============================================================
import sys as _sys, pathlib as _pathlib
_ROOT = next(_p for _p in _pathlib.Path(__file__).resolve().parents
             if (_p / 'config' / 'paths.py').exists())
_sys.path.insert(0, str(_ROOT))
import config  # noqa: F401  - also puts functions/ on sys.path
from config.paths import (
    ARCHIVE,
    PROCESSED_OUTPUT_DIR,
    ensure_dir,
)

Output_diri = str(ensure_dir(PROCESSED_OUTPUT_DIR)) + '/'
svante_archive = str(ARCHIVE) + '/'

# ============================================================

#================================================================================================
import os

MonitorInfo_filepath = f'{Output_diri}MonitorInfo/Lee_Jaffe_GAM_stats.csv'
Monitorne30Idx_filepath = f'{Output_diri}MonitorInfo/MatchedMonitors_ne30_ColIdx.csv'

# Variable Resolution Grid
SCRIP_ne30 = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../grid_files/ne30np4_091226_pentagons.nc')

# list for all case names
BASE2022_casename = 'f.e22.FCnudged.ne30_ne30_mg17.BGO3.BASEY20220401TY20230401'
BASE2023_casename = 'f.e22.FCnudged.ne30_ne30_mg17.BGO3.BASEY20230401TY20231101'
noBB2022_casename = 'f.e22.FCnudged.ne30_ne30_mg17.BGO3.noBBemisCONUS80kmBufferY20220401TY20221101'
noBB2023_casename = 'f.e22.FCnudged.ne30_ne30_mg17.BGO3.noBBemisCONUS80kmBufferY20230401TY20231101'
noAnthro2022_casename = 'f.e22.FCnudged.ne30_ne30_mg17.BGO3.noANTHROemisCONUS80kmBufferY20220401TY20221101'
noAnthro2023_casename = 'f.e22.FCnudged.ne30_ne30_mg17.BGO3.noANTHROemisCONUS80kmBufferY20230401TY20231101'


# label for each casename
casename_label_dic = {
    BASE2022_casename:'BASE2022',
    BASE2023_casename:'BASE2023',
    ### Perturbations
    noBB2022_casename:'noBB2022',
    noBB2023_casename:'noBB2023',
    noAnthro2022_casename:'noAnthro2022',
    noAnthro2023_casename:'noAnthro2023',
}

# merged h2 surface file for each case
mergedh2_surfO3filepath_dic = {
    BASE2022_casename:f'{Output_diri}h2_surfO3_merged/f.e22.FCnudged.ne30_ne30_mg17.BGO3.BASEY20220401TY20230401.cam.h2.surflev.O3.2022-04-01T2022-11-01.nc',
    BASE2023_casename:f'{Output_diri}h2_surfO3_merged/f.e22.FCnudged.ne30_ne30_mg17.BGO3.BASEY20230401TY20231101.cam.h2.surflev.O3.2023-04-01T2023-10-31.nc',
    ### Perturbations
    noBB2022_casename:f'{Output_diri}h2_surfO3_merged/f.e22.FCnudged.ne30_ne30_mg17.BGO3.noBBemisCONUS80kmBufferY20220401TY20221101.cam.h2.surflev.O3.2022-04-01T2022-10-31.nc',
    noBB2023_casename:f'{Output_diri}h2_surfO3_merged/f.e22.FCnudged.ne30_ne30_mg17.BGO3.noBBemisCONUS80kmBufferY20230401TY20231101.cam.h2.surflev.O3.2023-04-01T2023-11-01.nc',
    noAnthro2022_casename:f'{Output_diri}h2_surfO3_merged/f.e22.FCnudged.ne30_ne30_mg17.BGO3.noANTHROemisCONUS80kmBufferY20220401TY20221101.cam.h2.surflev.O3.2022-04-01T2022-11-01.nc',
    noAnthro2023_casename:f'{Output_diri}h2_surfO3_merged/f.e22.FCnudged.ne30_ne30_mg17.BGO3.noANTHROemisCONUS80kmBufferY20230401TY20231101.cam.h2.surflev.O3.2023-04-01T2023-11-01.nc',
}

#================================================================================================
### Specified input:
# For a given case
case_ls = [
    # BASE2022_casename,
    # BASE2023_casename,
    noBB2022_casename,
    noBB2023_casename,
    noAnthro2022_casename,
    noAnthro2023_casename,
]

startMMDD = '04-01'
endMMDD = '10-31'

#================================================================================================
# Import general functions
import os
import glob
import fnmatch

import pandas as pd
import geopandas as gpd
import xarray as xr
import numpy as np

import warnings
warnings.filterwarnings("ignore", category=FutureWarning, module="xarray.core.accessor_dt")

#================================================================================================
# Define functions specific for this application
# my functions
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../functions/'))
from SE_analysis import get_site_index

from timezonefinder import TimezoneFinder
from datetime import datetime
import pytz

"""Function to get the hour offset to UTC based on the location of the monitor"""
def get_summer_offset(lati: float, loni: float) -> int:
    # Find the timezone name from lat/lon
    tf = TimezoneFinder()
    tz_name = tf.timezone_at(lat=lati, lng=loni)
    if tz_name is None:
        raise ValueError("Could not determine timezone for given coordinates")

    # Pick a midsummer date (July 1) to ensure DST applies if relevant
    dt = datetime(datetime.now().year, 7, 1)

    # Localize to that timezone
    tz = pytz.timezone(tz_name)
    localized_dt = tz.localize(dt, is_dst=True)

    # Get offset from UTC in hours
    offset_hours = int(localized_dt.utcoffset().total_seconds() // 3600)
    return offset_hours
# # Example
# lati, loni = 33.545278, -86.549167
# print(get_summer_offset(lati, loni))  # → -5


"""Function to calculate MDA8 O3 in the local time following EPA's guidance"""
def compute_mda8_UTCoffset(
    o3_hourly_utc: xr.DataArray,
    datesv1: list[str],
    utc_offset_hours: int,         # e.g., -5 for CDT
    min_hours_per_block: int = 6,  # ≥6 of 8 hours to form a valid 8-h mean
    min_blocks_per_day: int = 13   # ≥13 of the 17 daily blocks required
) -> xr.DataArray:
    """
    Compute daily MDA8 O3 from hourly O3 timestamps in UTC,
    using a fixed integer local-time offset (hours) for the rolling windows.

    - Shifts the time coordinate by `utc_offset_hours` to a 'local' clock
      (local = UTC + offset), performs rolling and day grouping,
      then returns results indexed by local dates in `datesv1`.
    - Works with extra dims (e.g., 'site') and dask-backed arrays.
    """

    # 1) Shift time axis to 'local' clock (no tz objects, just a simple offset)
    t_local = pd.DatetimeIndex(o3_hourly_utc.time.values) + pd.to_timedelta(utc_offset_hours, "h")
    O3 = o3_hourly_utc.assign_coords(time=("time", t_local))
    
    # 2) 8-hour rolling mean with EPA min-hours criterion
    O3_8h = O3.rolling(time=8, min_periods=min_hours_per_block).mean()
    
    # 3) Keep only the 17 blocks ending 07..23 local
    O3_8h_17 = O3_8h.where(O3_8h["time"].dt.hour.isin(np.arange(7, 24)))

    # 4) Group by local civil day; compute daily max and valid-block count
    day_idx = pd.DatetimeIndex(O3_8h_17.time.values).floor("D")
    O3_8h_17 = O3_8h_17.assign_coords(day=("time", day_idx))
    
    mda8_by_day    = O3_8h_17.groupby("day").max("time", skipna=True)
    blocks_per_day = O3_8h_17.groupby("day").count("time")

    # 5) Apply day validity (≥13 of 17 blocks)
    valid = blocks_per_day >= min_blocks_per_day
    mda8_by_day = mda8_by_day.where(valid)

    # 6) Reindex to requested local dates and return with standard 'time' coord
    target_days = pd.to_datetime(datesv1)  # these are local dates
    out = mda8_by_day.reindex(day=target_days)
    out = out.rename({"day": "time"}).assign_coords(time=("time", target_days))
    return out
    
#================================================================================================
# Get the saved match column Idx
MonitorIdx_df = pd.read_csv(Monitorne30Idx_filepath)
# Remove rows where MUSICA0_colIndex is 'Find None' or NaN
MonitorIdx_df = MonitorIdx_df[MonitorIdx_df["MUSICA0_colIndex"] != "Find None"].copy()
MonitorIdx_df = MonitorIdx_df.dropna(subset=["MUSICA0_colIndex"])
# Now cast to int
MonitorIdx_df["MUSICA0_colIndex"] = MonitorIdx_df["MUSICA0_colIndex"].astype(int)

#================================================================================================
# Main function to calculate and save MDA8O3 and hourlyO3
def casei_build_and_save_all_sites(casename, MonitorIdx_df, startMMDD, endMMDD, Output_diri):
    """
    Build and save NetCDF datasets for daily local-time MDA8 O3 and hourly UTC O3
    at all monitor sites for one case/year, using AQS_code as the site dimension.
    """
    import pandas as pd
    import xarray as xr
    import numpy as np

    # ---- Build local-date range (inclusive) for target year ----
    YYYY = casename_label_dic[casename][-4:]
    startfileDate = f"{YYYY}-{startMMDD}"  # 'YYYY-MM-DD'
    endfileDate   = f"{YYYY}-{endMMDD}"    # 'YYYY-MM-DD'
    tdays = pd.date_range(start=startfileDate, end=endfileDate, freq="D")
    datesv1 = tdays.strftime("%Y-%m-%d").tolist()

    print(len(datesv1), "dates generated")
    print("First few v1:", datesv1[:5])

    # ---- Open merged hourly surface O3 (UTC timestamps) ----
    HourlyFilePath = mergedh2_surfO3filepath_dic[casename]
    HourlyO3_da = xr.open_dataset(HourlyFilePath)["O3"]  # expects dims incl. 'time' and 'ncol'

    # ---- Ensure AQS_code uniqueness and as string ----
    MonitorIdx_df = MonitorIdx_df.copy()
    MonitorIdx_df["AQS_code"] = MonitorIdx_df["AQS_code"].astype(str)
    dup_codes = MonitorIdx_df["AQS_code"][MonitorIdx_df["AQS_code"].duplicated(keep=False)]
    if not dup_codes.empty:
        raise ValueError(f"Duplicate AQS_code(s) found: {sorted(dup_codes.unique())}")

    # ---- Build per-site datasets with AQS_code as the DIMENSION ----
    mda8_list = []
    hourly_list = []

    for row in MonitorIdx_df.itertuples(index=False):
        colIdxi   = int(row.MUSICA0_colIndex)
        loni      = float(row.lon)
        lati      = float(row.lat)
        aqs_codei = str(row.AQS_code)

        # Hourly UTC series for this site (select by model column index)
        Monitori_HourlyO3_da = HourlyO3_da.sel(ncol=colIdxi)

        # Compute daily MDA8 for local dates using fixed summertime UTC offset
        MDA8O3 = compute_mda8_UTCoffset(
            Monitori_HourlyO3_da, datesv1,
            utc_offset_hours=get_summer_offset(lati, loni)
        )  # -> xr.DataArray with 'time' on local dates

        # Per-site MDA8 dataset (AQS_code x time)
        ds_mda8 = xr.Dataset(
            data_vars={
                "MDA8O3": (("AQS_code", "time"), MDA8O3.values[np.newaxis, :]*1e9, # to convert to ppb
                           {"units": "ppb",
                            "long_name": "Daily maximum of 8-hour average ozone (MDA8)",
                            "cell_methods": "time: mean (interval: 8 hours) time: maximum within days"})
            },
            coords={
                "AQS_code": ("AQS_code", [aqs_codei]),
                "time":     ("time", tdays),
                "lat":      ("AQS_code", [lati]),
                "lon":      ("AQS_code", [loni]),
            },
        )
        mda8_list.append(ds_mda8)

        # Per-site Hourly dataset (AQS_code x time) — UTC timestamps preserved
        ds_hourly = xr.Dataset(
            data_vars={
                "O3": (("AQS_code", "time"), Monitori_HourlyO3_da.values[np.newaxis, :],
                       {"units": Monitori_HourlyO3_da.attrs.get("units", "mol/mol"),
                        "long_name": "Hourly surface O3 (UTC)"})
            },
            coords={
                "AQS_code": ("AQS_code", [aqs_codei]),
                "time":     ("time", Monitori_HourlyO3_da.time.values),
                "lat":      ("AQS_code", [lati]),
                "lon":      ("AQS_code", [loni]),
            },
        )
        hourly_list.append(ds_hourly)

    # ---- Concatenate across monitors along AQS_code ----
    ds_all_mda8   = xr.concat(mda8_list,   dim="AQS_code").sortby("AQS_code")
    ds_all_hourly = xr.concat(hourly_list, dim="AQS_code").sortby("AQS_code")

    # ---- Coordinate/metadata niceties ----
    for ds in (ds_all_mda8, ds_all_hourly):
        ds["lat"].attrs.update({"units": "degrees_north", "standard_name": "latitude"})
        ds["lon"].attrs.update({"units": "degrees_east",  "standard_name": "longitude"})
        ds["AQS_code"].attrs.update({"long_name": "AQS site identifier"})

    ds_all_mda8.attrs.update({
        "title": "Daily MDA8 O3 at point monitors (local-time dates)",
        "Conventions": "CF-1.9",
        "comment": ("Local-time MDA8 computed from UTC hourly O3 using fixed summertime UTC offset "
                    "(8h rolling means with >=6 valid hours; day valid if >=13 of 17 end-hours 07–23)."),
    })
    ds_all_hourly.attrs.update({
        "title": "Hourly O3 at point monitors (UTC timestamps)",
        "Conventions": "CF-1.9",
        "comment": "Extracted from merged MUSICA outputs at monitor column indices (UTC timestamps).",
    })

    # ---- Save to NetCDF ----
    enc_mda8 = {
        "MDA8O3": {"zlib": True, "complevel": 5,
                   "chunksizes": (min(200, ds_all_mda8.dims["AQS_code"]), ds_all_mda8.dims["time"])},
        "lat": {"zlib": True, "complevel": 5},
        "lon": {"zlib": True, "complevel": 5},
        "AQS_code": {"zlib": True, "complevel": 5},
    }
    enc_hourly = {
        "O3": {"zlib": True, "complevel": 5,
               "chunksizes": (min(200, ds_all_hourly.dims["AQS_code"]), ds_all_hourly.dims["time"])},
        "lat": {"zlib": True, "complevel": 5},
        "lon": {"zlib": True, "complevel": 5},
        "AQS_code": {"zlib": True, "complevel": 5},
    }

    allMonitors_MDA8O3_filename = (
        f"{Output_diri}LocalTimeMDA8O3.{casename_label_dic[casename]}."
        f"GivenMonitors.{startfileDate}T{endfileDate}.nc"
    )
    allMonitors_HourlyO3_filename = (
        f"{Output_diri}UTChourlyO3.{casename_label_dic[casename]}."
        f"GivenMonitors.{startfileDate}T{endfileDate}.nc"
    )

    ds_all_mda8.to_netcdf(allMonitors_MDA8O3_filename, format="NETCDF4", encoding=enc_mda8)
    ds_all_hourly.to_netcdf(allMonitors_HourlyO3_filename, format="NETCDF4", encoding=enc_hourly)

    print(f"Saved MDA8 dataset to   {allMonitors_MDA8O3_filename}")
    print(f"Saved Hourly dataset to {allMonitors_HourlyO3_filename}")

#================================================================================================
# Apply this function to a selected list of case and provided dates
for casename in case_ls:
    print(f'Processing {casename}')
    #------------------------------
    saveto_diri = f'{Output_diri}ForGivenMonitors/'
    casei_build_and_save_all_sites(casename, MonitorIdx_df, startMMDD, endMMDD, saveto_diri)
print("Done!")