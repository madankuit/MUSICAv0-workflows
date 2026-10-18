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
- mda8_epa / lst_offset_hours (functions/mda8_epa.py):
    MDA8 per 40 CFR 50 Appendix U; local standard time offset of the monitor.

Outputs
-------
Two NetCDF files in Output_diri:
- LocalTimeMDA8O3.<case_label>.GivenMonitors.<startDate>T<endDate>.nc
    Daily MDA8 O3 at all monitor sites (units: ppb).
- UTChourlyO3.<case_label>.GivenMonitors.<startDate>T<endDate>.nc
    Hourly surface O3 at all monitor sites (units: model original).

Notes
-----
- MDA8 calculation: 40 CFR 50 Appendix U (functions/mda8_epa.py): local standard
  time, 8-h averages starting 07-23 (>=6 h, zero-substitution rule), truncation to
  whole ppb, day valid with >=13 of 17 averages or if the maximum exceeds 70 ppb.
- Output NetCDFs include per-site coordinates (lat, lon, AQS_code).
- Script prints the saved file paths; no return object.

MODIFICATION HISTORY:
    4 Sep 2025: VERSION 1.0
    - Initial version
    31 Aug 2026: VERSION 1.1
    - Paths, case names and the merged-file lookup moved to config/paths.py
    18 Oct 2026: VERSION 1.3
    - MDA8 exactly per 40 CFR 50 Appendix U via functions/mda8_epa.py (local
      standard time, windows starting 07-23, truncation, validity rules);
      --mda8-only, --cdate
    9 Oct 2026: VERSION 1.2
    - Loops over (scenario, year) and finds the merged file by (case, year):
      the noBBGlobal case is one run serving both years
    - --scenarios to process a subset; case/scenario/source file in attrs

    python Extract_givenmonitorO3_hourly_dailyMDA8_toNetCDF.py --scenarios noBBGlobal
    - --experiment ne0CONUSne30x8: same monitors and MDA8 code on the ne0CONUSne30x8 runs
      (column indices from GetMatched_..._ColumnIndex.py --experiment ne0CONUSne30x8)

    python Extract_givenmonitorO3_hourly_dailyMDA8_toNetCDF.py --experiment ne0CONUSne30x8
    - --match exact uses the exact point-in-polygon column indices (ne30: the
      _exact file); --file-tag goes after 'GivenMonitors' in output names so new
      files never overwrite delivered ones

    python Extract_givenmonitorO3_hourly_dailyMDA8_toNetCDF.py --match exact --file-tag All803
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
ap.add_argument('--match', choices=['legacy', 'exact'], default=None,
                help='column-index file: legacy (ne30 default, as delivered) or exact')
ap.add_argument('--file-tag', default='', help="inserted after 'GivenMonitors' in output names")
ap.add_argument('--mda8-only', action='store_true', help='write the MDA8 files only (hourly point files untouched)')
ap.add_argument('--cdate', default='', help="YYYYMMDD; MDA8 file names end in _c<cdate>.nc")
args = ap.parse_args()
EXP = BGO3_EXPERIMENTS[args.experiment]
CASES = EXP['cases']; LABELS = {k: f"{EXP['label_prefix']}{k[0]}{k[1]}" for k in CASES}
args.scenarios = args.scenarios or list(EXP['scenarios'])
assert set(args.scenarios) <= set(EXP['scenarios']), args.scenarios

MonitorInfo_filepath = BGO3_MONITOR_LIST
# column indices on this experiment's grid
_match = args.match or ('legacy' if args.experiment == 'ne30' else 'exact')
Monitorne30Idx_filepath = EXP['colidx_exact'] if (_match == 'exact' and 'colidx_exact' in EXP) else EXP['colidx']
print("column indices:", Monitorne30Idx_filepath)

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
from mda8_epa import mda8_epa, lst_offset_hours, REFERENCES as MDA8_REFERENCES
from timezonefinder import TimezoneFinder as _TFclass
_TF = _TFclass()

from timezonefinder import TimezoneFinder
from datetime import datetime
import pytz

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

        # Daily MDA8 per 40 CFR 50 Appendix U (local standard time, truncation, validity rules)
        _off = lst_offset_hours(_TF.timezone_at(lat=lati, lng=loni), loni, int(YYYY))
        _, _m = mda8_epa(Monitori_HourlyO3_da.values * 1e9, Monitori_HourlyO3_da.time.values, [_off],
                         startfileDate, endfileDate)
        mda8_ppb = _m[:, 0]                                   # whole ppb, NaN = not a valid day

        # Per-site MDA8 dataset (AQS_code x time)
        ds_mda8 = xr.Dataset(
            data_vars={
                "MDA8O3": (("AQS_code", "time"), mda8_ppb[np.newaxis, :],      # ppb
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
        "comment": "MDA8 per 40 CFR Part 50 Appendix U (2015 O3 NAAQS); see mda8_method.",
        "mda8_method": ("hourly values in local standard time (no daylight saving), truncated to whole ppb "
                        "(0.001 ppm); 8-h averages stored at their start hour, valid with >=6 of 8 hours (divisor "
                        "6/7), or with <=5 hours if the zero-substituted average exceeds 70 ppb; averages truncated "
                        "to whole ppb; daily maximum of the 17 averages starting 07:00-23:00 LST; day valid with "
                        ">=13 of 17 valid averages or if the maximum exceeds 70 ppb. Model hourly means are stamped "
                        "at the end of the hour, so hour-beginning = stamp - 1 h."),
        "mda8_values": "whole ppb (truncated, as EPA reports 3 decimals in ppm); NaN = not a valid day",
        "references": MDA8_REFERENCES,
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
        f"GivenMonitors{args.file_tag}.{startfileDate}T{endfileDate}"
        f"{'_c' + args.cdate if args.cdate else ''}.nc"
    )
    allMonitors_HourlyO3_filename = (
        f"{Output_diri}UTChourlyO3.{label}."
        f"GivenMonitors{args.file_tag}.{startfileDate}T{endfileDate}.nc"
    )

    assert not os.path.exists(allMonitors_MDA8O3_filename), allMonitors_MDA8O3_filename
    ds_all_mda8.to_netcdf(allMonitors_MDA8O3_filename, format="NETCDF4", encoding=enc_mda8)
    print(f"Saved MDA8 dataset to   {allMonitors_MDA8O3_filename}")
    if not args.mda8_only:
        ds_all_hourly.to_netcdf(allMonitors_HourlyO3_filename, format="NETCDF4", encoding=enc_hourly)
        print(f"Saved Hourly dataset to {allMonitors_HourlyO3_filename}")

#================================================================================================
# Apply this function to a selected list of case and provided dates
for scen, yr in run_ls:
    print(f'Processing {scen} {yr}')
    #------------------------------
    saveto_diri = str(ensure_dir(EXP['points_dir'])) + '/'
    casei_build_and_save_all_sites(scen, yr, MonitorIdx_df, startMMDD, endMMDD, saveto_diri)
print("Done!")