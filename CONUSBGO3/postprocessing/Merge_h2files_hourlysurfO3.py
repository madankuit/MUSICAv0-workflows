"""
    Merge the hourly h2 history files of every CONUSBGO3 case and keep the
    surface-layer O3 for a range of dates.

    All paths and case names come from config/paths.py; nothing is hard-coded.

MODIFICATION HISTORY:
    VERSION 1.0
    - Initial version
    31 Aug 2026: VERSION 1.1
    - Paths and case list moved to config/paths.py
    9 Oct 2026: VERSION 1.2
    - One merged file per (case, year), files picked by full date: the
      noBBGlobal case is a single run spanning both seasons
    - --scenarios to merge a subset; existing merged files are not redone
    - Checks: every file belongs to the case, no missing day, hourly time axis

    python Merge_h2files_hourlysurfO3.py --scenarios noBBGlobal
"""

# for hourly
histfreq  = "h2"

# specify the dates (MMDD, inclusive). Nov 1 is kept so local-time MDA8 on
# Oct 31 has its evening hours (local = UTC - 4..8 h).
startMMDD = "0401"
endMMDD   = "1101"

#================================================================================================
import numpy as np # for array manipulation and basic scientific calculation
import xarray as xr # To read NetCDF files
import matplotlib.pyplot as plt # Core library for plotting
import matplotlib.cm as cm # To use different colormaps
import cartopy.crs as ccrs # For map projection
import seaborn as sns # boxplot
import pandas as pd

# save PDF
from matplotlib.backends.backend_pdf import PdfPages

import warnings
# Ignore a specific warning
warnings.filterwarnings('ignore', message='Some warning message')
# Ignore the FutureWarning caused by iteritems()
warnings.filterwarnings("ignore", category=FutureWarning, message=".*iteritems.*")

import os
import re
import sys
import glob

### Box plot
import matplotlib.patches as mpatches ### , bbox_inches='tight'

#================================================================================================

# Configuration - every path and case name is imported from config/paths.py
import pathlib
_ROOT = next(p for p in pathlib.Path(__file__).resolve().parents
             if (p / 'config' / 'paths.py').exists())
sys.path.insert(0, str(_ROOT))
import config  # noqa: F401  - also puts functions/ on sys.path
from config.paths import (
    BGO3_CASES,
    BGO3_CASE_LABELS,
    BGO3_SCENARIOS,
    bgo3_merged_surfo3_glob,
    BGO3_MERGED_SURFO3_DIR,
    case_hist_dir,
    ensure_dir,
)

import argparse
ap = argparse.ArgumentParser()
ap.add_argument('--scenarios', nargs='+', default=list(BGO3_SCENARIOS),
                choices=list(BGO3_SCENARIOS))
args = ap.parse_args()

varlist = ['O3']
lev_idx = -1

ensure_dir(BGO3_MERGED_SURFO3_DIR)

#================================================================================================
### read in hourly mean
date_re = re.compile(r"\.(\d{4})-(\d{2})-(\d{2})-\d+\.nc$")  # ...YYYY-MM-DD-XXXXX.nc

def file_date(fname):
    m = date_re.search(fname)
    return pd.Timestamp(f"{m.group(1)}-{m.group(2)}-{m.group(3)}") if m else None

# Writing hourly data for the surface layer, one file per (case, year)
for (scen, yr), casename in BGO3_CASES.items():
    if scen not in args.scenarios:
        continue
    label = BGO3_CASE_LABELS[(scen, yr)]
    print(f"\n=== {label}: {casename}")

    t0 = pd.Timestamp(f"{yr}-{startMMDD[:2]}-{startMMDD[2:]}")
    t1 = pd.Timestamp(f"{yr}-{endMMDD[:2]}-{endMMDD[2:]}")
    RunPath = str(case_hist_dir(casename))
    all_files = glob.glob(os.path.join(RunPath, f"{casename}.cam.{histfreq}.*.nc"))
    RunFiles = sorted(f for f in all_files
                      if file_date(f) is not None and t0 <= file_date(f) <= t1)

    # ---- the selection must be this case, this year, every day ----
    days = pd.DatetimeIndex([file_date(f) for f in RunFiles])
    assert all(os.path.basename(f).startswith(casename + f".cam.{histfreq}.") for f in RunFiles)
    missing = pd.date_range(t0, min(t1, days.max())).difference(days)
    assert len(missing) == 0, f"{label}: missing h2 days {list(missing.strftime('%Y-%m-%d'))}"
    assert days.min() == t0 and days.max() >= pd.Timestamp(f"{yr}-10-31"), \
        f"{label}: files span {days.min().date()}..{days.max().date()}"
    print(f"Selected {len(RunFiles)} files {days.min().date()} .. {days.max().date()}")

    startfileDate = days.min().strftime('%Y-%m-%d')
    endfileDate = days.max().strftime('%Y-%m-%d')
    HourlyFilePath = str(BGO3_MERGED_SURFO3_DIR /
                         f'{casename}.cam.h2.surflev.O3.{startfileDate}T{endfileDate}.nc')
    if glob.glob(bgo3_merged_surfo3_glob(casename, yr)):
        print('Exists, not redone:', glob.glob(bgo3_merged_surfo3_glob(casename, yr)))
        continue

    # read in
    ds = xr.open_mfdataset(RunFiles, combine='nested', concat_dim=['time'], coords='minimal',
                           compat='override', use_cftime=False)
    dt = np.diff(ds.time.values).astype('timedelta64[m]').astype(int)
    assert (dt == 60).all(), f"{label}: time axis not hourly (steps {sorted(set(dt))})"
    # surface
    surf_da = ds.isel(lev=lev_idx, ilev=lev_idx)['O3']
    surf_da.attrs.update(case_name=casename, scenario=scen, year=yr,
                         source_files=f"{os.path.basename(RunFiles[0])} .. {os.path.basename(RunFiles[-1])}")

    # save to .nc
    surf_da.to_netcdf(HourlyFilePath)
    print('Save to:', HourlyFilePath)
