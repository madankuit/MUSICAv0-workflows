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
sys.path.insert(0,'/home/taoma528/Scripts/CESM_analysis/functions/')

from func_MUSICA_DefineRegion import *
from func_ModelEval_statistical_tests import *
from func_MUSICAvsSLAMS_surf import *
# from func_ReadMUSICAOutput import *

from Plot_2D import Plot_2D # To draw a map
from SE_analysis import get_site_index

#================================================================================================

svante_archive = '/net/fs09/d0/taoma528/CESM22/archive/'
varlist = ['area','lat','lon','CH2O','NO2','O3']
lev_idx = -1

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
OnlyNOHourlyNEI_casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1basehourlyNOotherJulyMeanNEI2017'
HourlyNEI_casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1basehourlyNEI2017'
nudge6HourlyNEI_casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.6HrNudgeTS1hourlyNEI2017'

casename_ls = [
                Base_casename,
               # MonthlyNEI_casename,
               # OnlyNOHourlyNEI_casename,
               # HourlyNEI_casename,
               # nudge6HourlyNEI_casename,
              ]

# label for each casename
casename_label_dic = {
                'SLAMS':'SLAMS',
                'TROPOMI':'TROPOMI',
                Base_casename:'Base',
                MonthlyNEI_casename:'Monthly NEI',
                OnlyNOHourlyNEI_casename:'Hourly-NO NEI',
                HourlyNEI_casename:'Hourly NEI',
                nudge6HourlyNEI_casename:'6-hr Nudge - Hourly NEI',
                }

# Define colors for each case
casename_color_dic = {
                'SLAMS':'black',
                'TROPOMI':'black',
                Base_casename:'tab:blue',
                MonthlyNEI_casename:'tab:olive',
                OnlyNOHourlyNEI_casename:'salmon',
                HourlyNEI_casename:'tab:purple',
                nudge6HourlyNEI_casename:'tab:cyan',
                     }
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
    HourlyFilePath = f'/net/fs09/d0/taoma528/CESM22/processed_output/July2018_surfh2_merged/{casename}.cam.h2.surflev.{startfileDate}T{endfileDate}.nc'
    
    # save to .nc
    surf_ds.to_netcdf(HourlyFilePath)
    print('Save to:',HourlyFilePath)

#================================================================================================


#================================================================================================
