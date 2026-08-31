"""
    Merge the hourly h2 history files of every CONUSBGO3 case and keep the
    surface-layer O3 for a range of dates.

    All paths and case names come from config/paths.py; nothing is hard-coded.

MODIFICATION HISTORY:
    VERSION 1.0
    - Initial version
    31 Aug 2026: VERSION 1.1
    - Paths and case list moved to config/paths.py
"""

# for hourly
histfreq  = "h2"

# specify the dates (MMDD, inclusive)
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
    BGO3_MERGED_SURFO3_DIR,
    case_hist_dir,
    ensure_dir,
)

varlist = ['O3']
lev_idx = -1

casename_ls = list(BGO3_CASES.values())
casename_label_dic = dict(BGO3_CASE_LABELS)

ensure_dir(BGO3_MERGED_SURFO3_DIR)
                     
#================================================================================================
### read in hourly mean
histfreq = 'h2'

# Writing hourly data for the surface layer 
for caseIdx in range(len(casename_ls)):
    casename = casename_ls[caseIdx]
    print(casename)
    # year information is contained in the label name as the last four digits
    # YYYY = casename_label_dic[casename][-4:]
    
    # select the files
    RunPath = str(case_hist_dir(casename))
    date_re = re.compile(r"\.(\d{4})-(\d{2})-(\d{2})-\d+\.nc$")  # ...YYYY-MM-DD-XXXXX.nc

    def extract_mmdd(fname: str) -> str:
        m = date_re.search(fname)
        if not m:
            return None
        # Return 'MMDD'
        return f"{m.group(2)}{m.group(3)}"

    def mmdd_in_range(mmdd: str, start_mmdd: str, end_mmdd: str) -> bool:
        """Inclusive range check on MMDD; supports wrap-around (e.g., Nov–Mar)."""
        if mmdd is None:
            return False
        if start_mmdd <= end_mmdd:
            return start_mmdd <= mmdd <= end_mmdd
        else:
            # wrap-around across year end (e.g., start='1101', end='0201')
            return (mmdd >= start_mmdd) or (mmdd <= end_mmdd)

    # Gather and filter
    all_files = glob.glob(os.path.join(RunPath, f"*{histfreq}*.nc"))
    RunFiles  = sorted(f for f in all_files if mmdd_in_range(extract_mmdd(f), startMMDD, endMMDD))
    print(f"Selected {len(RunFiles)} files in MMDD [{startMMDD}–{endMMDD}]")

    
    # read in 
    ds = xr.open_mfdataset(RunFiles, combine='nested', concat_dim=['time'], coords='minimal', compat='override', use_cftime=False)
    # surface
    surf_da = ds.isel(lev=lev_idx,ilev=lev_idx)['O3']
    
    startfileDate = RunFiles[0].split('.')[-2][:10]
    endfileDate = RunFiles[-1].split('.')[-2][:10]
    # Save
    HourlyFilePath = str(BGO3_MERGED_SURFO3_DIR /
                         f'{casename}.cam.h2.surflev.O3.{startfileDate}T{endfileDate}.nc')
    
    # save to .nc
    surf_da.to_netcdf(HourlyFilePath)
    print('Save to:',HourlyFilePath)

#================================================================================================


#================================================================================================
