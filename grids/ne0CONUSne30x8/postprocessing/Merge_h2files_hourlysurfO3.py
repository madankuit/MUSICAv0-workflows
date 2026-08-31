"""
    This script merges the h2 files to get the surface O3 for a range of dates
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

svante_archive = str(ARCHIVE) + '/'
Output_diri = str(ensure_dir(PROCESSED_OUTPUT_DIR)) + '/'

# ============================================================

# for hourly
histfreq  = "h2"

# specify the dates
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
import glob

### Box plot
import matplotlib.patches as mpatches ### , bbox_inches='tight'

#================================================================================================

varlist = ['O3']
lev_idx = -1

# list for all case names
BASE2022_casename = 'f.e22.FCnudged.ne30_ne30_mg17.BGO3.BASEY20220401TY20230401'
BASE2023_casename = 'f.e22.FCnudged.ne30_ne30_mg17.BGO3.BASEY20230401TY20231101'
noBB2022_casename = 'f.e22.FCnudged.ne30_ne30_mg17.BGO3.noBBemisCONUS80kmBufferY20220401TY20221101'
noBB2023_casename = 'f.e22.FCnudged.ne30_ne30_mg17.BGO3.noBBemisCONUS80kmBufferY20230401TY20231101'
noAnthro2022_casename = 'f.e22.FCnudged.ne30_ne30_mg17.BGO3.noANTHROemisCONUS80kmBufferY20220401TY20221101'
noAnthro2023_casename = 'f.e22.FCnudged.ne30_ne30_mg17.BGO3.noANTHROemisCONUS80kmBufferY20230401TY20231101'

casename_ls = [
               BASE2022_casename,
               BASE2023_casename,
               noBB2022_casename,
               noBB2023_casename,
               noAnthro2022_casename,
               noAnthro2023_casename,
              ]

# label for each casename
casename_label_dic = {
                'SLAMS':'SLAMS',
                BASE2022_casename:'BASE2022',
                BASE2023_casename:'BASE2023',
                ### Perturbations
                noBB2022_casename:'noBB2022',
                noBB2023_casename:'noBB2023',
                noAnthro2022_casename:'noAnthro2022',
                noAnthro2023_casename:'noAnthro2023',
                }
                     
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
    # RunPath = f'{svante_archive}{casename}/atm/hist/'
    # RunFiles = sorted(glob.glob(os.path.join(RunPath, f'*{histfreq}*')))
    
    RunPath = f"{svante_archive}{casename}/atm/hist/"
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
    HourlyFilePath = f'{Output_diri}h2_surfO3_merged/{casename}.cam.h2.surflev.O3.{startfileDate}T{endfileDate}.nc'
    
    # save to .nc
    surf_da.to_netcdf(HourlyFilePath)
    print('Save to:',HourlyFilePath)

#================================================================================================


#================================================================================================
