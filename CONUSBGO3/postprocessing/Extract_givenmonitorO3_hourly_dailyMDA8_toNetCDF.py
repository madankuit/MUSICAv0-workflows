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
    Directory where output NetCDF files are written (config.paths
    BGO3_GIVEN_MONITORS_DIR).
- BGO3_CASE_LABELS: dict
    Maps (scenario, year) to string label, e.g. 'noBBGlobal2023'.
- bgo3_merged_surfo3_path(casename, year):
    Merged hourly surface O3 file for one case and year.
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
    4 Sep 2025: VERSION 1.0
    - Initial version
    31 Aug 2026: VERSION 1.1
    - Paths, case names and the merged-file lookup moved to config/paths.py
    9 Oct 2026: VERSION 1.2
    - Loops over (scenario, year) and finds the merged file by (case, year):
      the noBBGlobal case is one run serving both years
    - --scenarios to process a subset; case/scenario/source file in attrs

    python Extract_givenmonitorO3_hourly_dailyMDA8_toNetCDF.py --scenarios noBBGlobal
    - --experiment ne0CONUS: same monitors and MDA8 code on the ne0CONUSne30x8 runs
      (column indices from GetMatched_..._ColumnIndex.py --experiment ne0CONUS)

    python Extract_givenmonitorO3_hourly_dailyMDA8_toNetCDF.py --experiment ne0CONUS
"""
#================================================================================================
#================================================================================================
# Configuration - every path and case name is imported from config/paths.py
import sys
import pathlib
_ROOT = next(p for p in pathlib.Path(__file__).resolve().parents
             if (p / 'config' / 'paths.py').exists())
sys.path.insert(0, str(_ROOT))
import config  # noqa: F401  - also puts functions/ on sys.path
from config.paths import (
    BGO3_CASES,
    BGO3_CASE_LABELS,
    BGO3_SCENARIOS,
    BGO3_EXPERIMENTS,
    BGO3_MONITOR_LIST,
    BGO3_MONITOR_COLIDX,
    BGO3_GIVEN_MONITORS_DIR,
    BGO3_START_MMDD,
    BGO3_END_MMDD,
    SCRIP_NE30NP4,
    bgo3_merged_surfo3_path,
    ensure_dir,
)

import argparse
ap = argparse.ArgumentParser()
ap.add_argument('--experiment', default='ne30', choices=list(BGO3_EXPERIMENTS))
ap.add_argument('--scenarios', nargs='+', default=None)
args = ap.parse_args()
EXP = BGO3_EXPERIMENTS[args.experiment]
CASES = EXP['cases']; LABELS = {k: f"{EXP['label_prefix']}{k[0]}{k[1]}" for k in CASES}
args.scenarios = args.scenarios or list(EXP['scenarios'])
assert set(args.scenarios) <= set(EXP['scenarios']), args.scenarios

MonitorInfo_filepath = BGO3_MONITOR_LIST
Monitorne30Idx_filepath = EXP['colidx']      # column indices on this experiment's grid

# Variable Resolution Grid (ships with the repo)
SCRIP_ne30 = str(SCRIP_NE30NP4)

#================================================================================================
### Specified input:
# (scenario, year) pairs to process; default is every CONUSBGO3 case.
# The merged file is resolved per (case, year) by bgo3_merged_surfo3_path.
run_ls = [k for k in CASES if k[0] in args.scenarios]

startMMDD = BGO3_START_MMDD
endMMDD = BGO3_END_MMDD

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
# repo-local shared functions (functions/ is on sys.path via `import config`)
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
def casei_build_and_save_all_sites(scen, yr, MonitorIdx_df, startMMDD, endMMDD, Output_diri):
    """
    Build and save NetCDF datasets for daily local-time MDA8 O3 and hourly UTC O3
    at all monitor sites for one (scenario, year), using AQS_code as the site dimension.
    """
    import pandas as pd
    import xarray as xr
    import numpy as np

    # ---- Build local-date range (inclusive) for target year ----
    casename = CASES[(scen, yr)]
    label = LABELS[(scen, yr)]
    YYYY = str(yr)
    startfileDate = f"{YYYY}-{startMMDD}"  # 'YYYY-MM-DD'
    endfileDate   = f"{YYYY}-{endMMDD}"    # 'YYYY-MM-DD'
    tdays = pd.date_range(start=startfileDate, end=endfileDate, freq="D")
    datesv1 = tdays.strftime("%Y-%m-%d").tolist()

    print(len(datesv1), "dates generated")
    print("First few v1:", datesv1[:5])

    # ---- Open merged hourly surface O3 (UTC timestamps) ----
    HourlyFilePath = bgo3_merged_surfo3_path(casename, yr, EXP['merged_dir'])
    print(f"{label}: {casename}\n  merged file: {HourlyFilePath}")
    HourlyO3_da = xr.open_dataset(HourlyFilePath)["O3"]  # expects dims incl. 'time' and 'ncol'
    assert HourlyO3_da.sizes["ncol"] == EXP["ncol"], (HourlyFilePath, HourlyO3_da.sizes["ncol"])
    _t = pd.DatetimeIndex(HourlyO3_da.time.values)
    assert (_t.year == yr).all() and _t.min() <= pd.Timestamp(f"{startfileDate} 01:00") \
        and _t.max() >= pd.Timestamp(f"{endfileDate} 23:00"), \
        f"{label}: merged file spans {_t.min()} .. {_t.max()}"

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

    _src = {"case_name": casename, "scenario": scen, "year": yr, "grid": EXP["grid"],
            "model": EXP["model"], "column_index_file": os.path.basename(str(Monitorne30Idx_filepath)),
            "merged_input_file": os.path.basename(HourlyFilePath)}
    ds_all_mda8.attrs.update(_src)
    ds_all_hourly.attrs.update(_src)
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
        f"{Output_diri}LocalTimeMDA8O3.{label}."
        f"GivenMonitors.{startfileDate}T{endfileDate}.nc"
    )
    allMonitors_HourlyO3_filename = (
        f"{Output_diri}UTChourlyO3.{label}."
        f"GivenMonitors.{startfileDate}T{endfileDate}.nc"
    )

    ds_all_mda8.to_netcdf(allMonitors_MDA8O3_filename, format="NETCDF4", encoding=enc_mda8)
    ds_all_hourly.to_netcdf(allMonitors_HourlyO3_filename, format="NETCDF4", encoding=enc_hourly)

    print(f"Saved MDA8 dataset to   {allMonitors_MDA8O3_filename}")
    print(f"Saved Hourly dataset to {allMonitors_HourlyO3_filename}")

#================================================================================================
# Apply this function to a selected list of case and provided dates
for scen, yr in run_ls:
    print(f'Processing {scen} {yr}')
    #------------------------------
    saveto_diri = str(ensure_dir(EXP['points_dir'])) + '/'
    casei_build_and_save_all_sites(scen, yr, MonitorIdx_df, startMMDD, endMMDD, saveto_diri)
print("Done!")