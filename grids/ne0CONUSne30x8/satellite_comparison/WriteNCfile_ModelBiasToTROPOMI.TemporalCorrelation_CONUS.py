# This script is used for CESM2.2 MUSICA simulations & TROPOMI proceesed to 0.15x0.15 degree resolution
'''
This code is designed to calculate Spearman correlation, Normalized Model Bias and Normalized RMSE for vertical column densities of HCHO, NO2, and CO of h2 MUSICA-V0 outputs approximated at 1:30 PM LT in a regular 0.15x0.15

MODIFICATION HISTORY:
    Match, 23, 2024: VERSION 1.0
    - Initial version
    
'''
#================================================================================================

varlist = ['NO2','HCHO','CO']

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
    CESM_FIGURE_DIR,
    REGRIDDED_2018_ROOT,
    TROPOMI_015_DIRS,
    ensure_dir,
)

svante_archive = str(ARCHIVE) + '/'
figure_diri = str(ensure_dir(CESM_FIGURE_DIR / 'ModelBiasToTROPOMI')) + '/'
L3ProductFileOut_diri_dic = {k: str(v) + '/' for k, v in TROPOMI_015_DIRS.items()}
MUSICA_regridded_base_diri = str(REGRIDDED_2018_ROOT) + '/'

# ============================================================

### Designated path and files
#-----------------------------------------------------------------------
# MUSICA
import os
SCRIP_CONUS = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../grid_files/ne0CONUS_ne30x8_np4_SCRIP.nc')

fileregion_ls = ['WestCoast','Mountain','Midwest','Southwest','Northeast','Southeast']
nfileregion = len(fileregion_ls)

Base_casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1base01'
MonthlyNEI_casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1baseJulyMeanNEI2017'
OnlyNOHourlyNEI_casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1basehourlyNOotherJulyMeanNEI2017'
HourlyNEI_casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1basehourlyNEI2017'
nudge6HourlyNEI_casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.6HrNudgeTS1hourlyNEI2017'

casename_ls = [
               Base_casename,
               MonthlyNEI_casename,
               OnlyNOHourlyNEI_casename,
               HourlyNEI_casename,
               # nudge6HourlyNEI_casename,
              ]

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

# label for each casename
casename_label_dic = {
                'SLAMS':'SLAMS',
                'TROPOMI':'TROPOMI',
                Base_casename:'Base',
                MonthlyNEI_casename:'Monthly NEI',
                OnlyNOHourlyNEI_casename:'Hourly NO NEI',
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

# Get TROPOMI data
fileheader_dic = {'HCHO':'S5P_RPRO_L2__HCHO___',
                   'NO2':'S5P_RPRO_L2__NO2____',
                   'CO':'S5P_RPRO_L2__CO_____',
                  }

# L3ProductFileOut_diri_dic is set in the USER CONFIGURATION block above

#================================================================================================
#### Module import ###
import os
import glob

import numpy as np # for array manipulation and basic scientific calculation
import xarray as xr # To read NetCDF files
import pandas as pd

import matplotlib.pyplot as plt

import warnings
# Ignore a specific warning
warnings.filterwarnings('ignore', message='Some warning message')
# Ignore the FutureWarning caused by iteritems()
warnings.filterwarnings("ignore", category=FutureWarning, message=".*iteritems.*")

#================================================================================================
### Self-defined functions
# Calculate the mean concentrations of each species averaged at each region
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../functions/'))
from func_MUSICA_DefineRegion import *
from func_ModelEval_statistical_tests import *
# SCRIP_CONUS and figure_diri are set in the USER CONFIGURATION block above

#================================================================================================
### Date for July 2018
startdate = "2018-07-01" #"2018-06-30"
enddate = "2018-08-01" #"2018-08-03"
filetime = f'{startdate}T{enddate}'

# Date format dictionary for processing dates
# Find the list of dates
from datetime import datetime, timedelta
from pytz import timezone

TROPOMI_days = []
MUSICA_days = []
import datetime
Read_start = datetime.datetime.strptime(startdate, "%Y-%m-%d")
Read_end_plus1 = datetime.datetime.strptime(enddate, "%Y-%m-%d")
date_generated = [Read_start + timedelta(days=x) for x in range(0, (Read_end_plus1-Read_start).days)]

# In the format of filename date
for date in date_generated:
    TROPOMI_days.append(date.strftime("%Y%m%d"))
    MUSICA_days.append(date.strftime("%Y-%m-%d"))
    
# in a dictionary to adjust for different date format
TRROPOMIdate_onMUSICAdate = dict(zip(TROPOMI_days,MUSICA_days))

#================================================================================================
#================================================================================================
#### For Daily (MUSICA uses daily mean, while TROPOMI is one day)
#### Calculate r, NMBE, NRMSE at each pixel, save to .nc file 
#-----------------------------------------------------------------------
### FOR EACH VAR
#-----------------------------------------------------------------------        
for varname in varlist:
    #-----------------------------------------------------------------------
    ### FOR EACH CASE: Get the values for each region
    #-----------------------------------------------------------------------
    for casename in casename_ls:
        
        if varname in ['NO2','HCHO']:
            colOpt = 'TropVCD'
            colLabel = r'$VCD_{Trop}$'
            FileOut_diri = f'{MUSICA_regridded_base_diri}{casename}/TropVCD_approx1330LT/'
        elif varname in ['CO']:
            colOpt = 'TotalVCD'
            colLabel = r'$VCD_{Total}$'
            FileOut_diri = f'{MUSICA_regridded_base_diri}{casename}/TotalVCD_approx1330LT/'

        #Save for each case: where to store the data
        SavePath = f'{FileOut_diri}{casename}.ModelBiasToTROPOMI.TemporalCor.CONUS.nc'
        #-----------------------------------------------------
        if os.path.exists(SavePath):
            print(f"The file {SavePath} exists.")
        else:
            print(f'Start to process for {varname}, for {casename}')
            ### Proceed
            ### Get TROPOMI for the given variable
            # Read in a list of TROPOMI filenames
            TROPOMI_VCD_filename_ls = []
            for datei in TROPOMI_days:
                filepath = L3ProductFileOut_diri_dic[varname]+fileheader_dic[varname]+'Global_Regrid015deg_'+datei+'.nc'
                if os.path.exists(filepath):
                    TROPOMI_VCD_filename_ls.append(filepath)
                else:
                    print(f'Missing TROPOMI file on {datei}')        
            # Global: Open all TROPOMI files and create a dataset
            JulyTROPOMI_VCD_ds = xr.open_mfdataset(TROPOMI_VCD_filename_ls)
            Multifactor = JulyTROPOMI_VCD_ds['TROPOMI_'+varname+'_'+colOpt].multiplication_factor_to_convert_to_molecules_percm2
            JulyDaily_TROPOMI_VCD_ds = JulyTROPOMI_VCD_ds['TROPOMI_'+varname+'_'+colOpt]
            #-----------------------------------------------------
            for fileregionidx in range(len(fileregion_ls)):
                fileregion = fileregion_ls[fileregionidx]
                # print(fileregion)
                lon_right,lat_bot,lon_left,lat_up = latlonbound(fileregion)
                #-----------------------------------------------------------------------------------
                # Select TROPOMI for the given region
                rawregioni_JulyDaily_TROPOMI_VCD_ds = JulyDaily_TROPOMI_VCD_ds.sel(lat=slice(lat_bot,lat_up),lon=slice(lon_right,lon_left))
                regioni_JulyDaily_TROPOMI_VCD_ds = rawregioni_JulyDaily_TROPOMI_VCD_ds*Multifactor
                #-----------------------------------------------------------------------------------
                # Select MUSICA for the given region | use maksed
                MUSICAVCDSave_diri = f'{MUSICA_regridded_base_diri}{casename}/{colOpt}_approx1330LT/'
                regioni_MUSICA_filename = f'{MUSICAVCDSave_diri}MaskedTROPOMInan_{casename}.{colOpt}1330LT_withTROPOMIAK.{varname}.{fileregion}.{filetime}.nc'
                regioni_MUSICA_VCD_ds = xr.open_dataset(regioni_MUSICA_filename)
                # MUSICA July daily mean 
                regioni_JulyDaily_MUSICA_VCD_ds = regioni_MUSICA_VCD_ds[varname+'_'+colOpt]

                ## dimensions
                time_ar = regioni_JulyDaily_TROPOMI_VCD_ds.time.values
                lat_ar = regioni_JulyDaily_TROPOMI_VCD_ds.lat.values
                lon_ar = regioni_JulyDaily_TROPOMI_VCD_ds.lon.values

                ## prepare array to store the correlation factors
                r_ar = np.zeros([len(lat_ar),len(lon_ar)])
                r_ar[:,:] = np.nan
                NMBE_ar = np.zeros([len(lat_ar),len(lon_ar)])
                NMBE_ar[:,:] = np.nan
                NRMSE_ar = np.zeros([len(lat_ar),len(lon_ar)])
                NRMSE_ar[:,:] = np.nan

                ## For each location (lat,lon)
                for latIdx in range(len(lat_ar)):
                    lati  = lat_ar[latIdx]
                    for lonIdx in range(len(lon_ar)):
                        loni  = lon_ar[lonIdx]
                        ### Get the Model and Observation values
                        MUSICAvals = regioni_JulyDaily_MUSICA_VCD_ds.isel(lat=latIdx,lon=lonIdx).values
                        TROPOMIvals = regioni_JulyDaily_TROPOMI_VCD_ds.isel(lat=latIdx,lon=lonIdx).values
                        ### correlation
                        r_ar[latIdx,lonIdx] = spearmanr_TwoArComp_returnr(MUSICAvals, TROPOMIvals)
                        NMBE_ar[latIdx,lonIdx] = calculate_NMBE(MUSICAvals, TROPOMIvals)
                        NRMSE_ar[latIdx,lonIdx] = calculate_NRMSE(MUSICAvals, TROPOMIvals)

                # write this variable to a dataarray
                casei_ds = xr.Dataset(
                                    {
                                        'spearmanr': (['lat', 'lon'], r_ar),
                                        'NMBE': (['lat', 'lon'], NMBE_ar),
                                        'NRMSE': (['lat', 'lon'], NRMSE_ar),
                                    },
                                    coords={'lat': ('lat', lat_ar), 
                                            'lon': ('lon', lon_ar)},
                                    attrs={'description': 'Model evaluation for daily July',
                                          'region': fileregion,},
                                    )

                # Add the attributes of coordinates
                for var_name in ['lat','lon']:    
                    # Copy attributes from source to target variable
                    for attr_name, attr_value in regioni_JulyDaily_MUSICA_VCD_ds[var_name].attrs.items():
                        casei_ds[var_name].attrs[attr_name] = attr_value

                # Merge for all regions 
                if fileregionidx==0:
                    merged_ds = casei_ds
                else:
                    # Combine the datasets for different regions
                    merged_ds = merged_ds.combine_first(casei_ds)

            #================================================================================================
            # save to .nc
            merged_ds.to_netcdf(SavePath)
            print('Save to:',SavePath)
