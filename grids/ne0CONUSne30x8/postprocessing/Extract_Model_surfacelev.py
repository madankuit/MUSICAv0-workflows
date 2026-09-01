### Extract_Model_surfacelev
### This script extracts the surface level

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
    JULY2018_SURFH2_MERGED_DIR,
    ensure_dir,
)

svante_archive = str(ARCHIVE) + '/'
processed_output_diri = str(ensure_dir(JULY2018_SURFH2_MERGED_DIR)) + '/'

# ============================================================

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
import glob

### Box plot
import matplotlib.patches as mpatches ### , bbox_inches='tight'

#================================================================================================
# functions I defined
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../functions/'))

from func_MUSICA_DefineRegion import *
# from func_ReadMUSICAOutput import *

from Plot_2D import Plot_2D # To draw a map
from SE_analysis import get_site_index

#================================================================================================

lev_idx = -1 # for near-surface level

varlabel_dic = {'O3':r'$O_{3}$',
                'NO2':r'$NO_{2}$',
                'HCHO':r'$HCHO$', # listing twice if pointed
               'CH2O':r'$HCHO$', #r'$CH_{2}O$'
                'CO':r'$CO$',
                'SO2':r'$SO_{2}$',
                'PM25':r'$PM_{2.5}$',
               }

unit_dic = {'O3':r'$ppb$',
                'NO2':r'$ppb$',
            'HCHO':r'$ppb$', # listing twice if pointed
               'CH2O':r'$ppb$',
                'CO':r'$ppb$',
                'SO2':r'$ppb$',
                'PM25':r'$\mu g/m^{3}$',   #r'$kg/m^{3}$', 
               }

Scalefactor_dic = {'O3':1e9,
                'NO2':1e9,
               'CH2O':1e9,
                'CO':1e9,
                'SO2':1e9,
                'PM25':1e9,
               }

fileregion_ls = ['WestCoast','Mountain','Midwest','Southwest','Northeast','Southeast']
nfileregion = len(fileregion_ls)

Base_casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1base01'
MonthlyNEI_casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1baseJulyMeanNEI2017'

HourlyNEI_casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1basehourlyNEI2017'

OnlyNOHourlyNEI_casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1basehourlyNOotherJulyMeanNEI2017'
OnlyNODailyNEI_casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1basedailyNOotherJulyMeanNEI2017'

nudge6HourlyNEI_casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.6HrNudgeTS1hourlyNEI2017'

###
v2MonthlyNEI_casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1baseJulyMeanNEI2017v2'
# perturbation cases
MonthlyNEIm30anthroNO_casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1baseJulyMeanNEI2017_m30anthroNO'
MonthlyNEIm30anthroVOC_casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1baseJulyMeanNEI2017_m30anthroVOC'


casename_ls = [
    # Base_casename,
   # MonthlyNEI_casename,
   # HourlyNEI_casename,
   # OnlyNOHourlyNEI_casename,
   # OnlyNODailyNEI_casename,
   # nudge6HourlyNEI_casename,
    # v2MonthlyNEI_casename,
    # MonthlyNEIm30anthroNO_casename,
    MonthlyNEIm30anthroVOC_casename,
              ]

#================================================================================================
### read in hourly mean
histfreq = 'h2'

# Writing hourly data for the surface layer 
for caseIdx in range(len(casename_ls)):
    casename = casename_ls[caseIdx]
    print(casename)
    
    # read in 
    RunPath = f'{svante_archive}{casename}/atm/hist/'
    RunFiles = sorted(glob.glob(os.path.join(RunPath, f'*{histfreq}*')))
    ds = xr.open_mfdataset(RunFiles, combine='nested', concat_dim=['time'], coords='minimal', compat='override', use_cftime=False)
    # surface
    surf_ds = ds.isel(lev=lev_idx,ilev=lev_idx)
    
    startfileDate = RunFiles[0].split('.')[-2][:10]
    endfileDate = RunFiles[-1].split('.')[-2][:10]
    # Save
    HourlyFilePath = f'{processed_output_diri}{casename}.cam.h2.surflev.{startfileDate}T{endfileDate}.nc'
    
    # save to .nc
    surf_ds.to_netcdf(HourlyFilePath)
    print('Save to:',HourlyFilePath)

#================================================================================================


#================================================================================================
