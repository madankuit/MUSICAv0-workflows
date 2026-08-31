# This script is used for CESM2.2 CAM-chem with SE dynamical core simulations & TROPOMI proceesed to 0.15x0.15 degree resolution
'''
this code is designed to calculate tropospheric vertical column densities of HCHO or NO2 of h2 MUSICA-V0 outputs approximated at 1:30 PM LT in a regular 0.15x0.15 lat and lon grids using the averaging kernels from TROPOMI

process separately for HCHO and NO2, for the duration of the given period

MODIFICATION HISTORY:
    December, 13, 2023: VERSION 1.0
    - Initial version
    December, 16, 2023: VERSION 1.1
    - Modify the way of reading in files to minimize memory demand,, adjust for different
    Feb, 25, 2024: VERSION 1.0
    - Bug fix, add one to the last layer index as slice will not include the right-hand number: indices[-1]+1
'''
#================================================================================================
# inputs
varname = 'HCHO' # 'HCHO' or 'NO2'

fileregion_ls = ['WestCoast','Mountain','Midwest','Southwest','Southeast','Northeast']
#Options: 'WestCoast','Mountain','Midwest','Southwest','Southeast','Northeast'

# Variables to modify
startdate = "2018-07-01" #"2018-06-30"
enddate = "2018-08-01" #"2018-08-01"

# casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1base01'
# casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1basehourlyNEI2017'
# casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1basehourlyNOotherJulyMeanNEI2017'
# casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1baseJulyMeanNEI2017'
### with 6-hr nudging
casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.6HrNudgeTS1hourlyNEI2017'

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
    REGRID_SUBDIR_BYVAR,
    REGRID_SUBDIR_TROPVCD,
    case_regrid_dir,
    ensure_dir,
)

# Both depend on `casename`, which each script sets below.
regridByVar_diri = ''
MUSICATropVCDSave_diri = ''
# Resolve them once `casename` is known:
#   regridByVar_diri = str(case_regrid_dir(casename, REGRID_SUBDIR_BYVAR)) + '/'
#   MUSICATropVCDSave_diri = str(ensure_dir(case_regrid_dir(casename, REGRID_SUBDIR_TROPVCD))) + '/'

# ============================================================

print(regridByVar_diri)

### Names of the files
TROP_P_filename = casename+'.cam.h2.MassConserve_latlon015.TROP_P.20180701T20180801.nc'
PMID_filename = casename+'.cam.h2.MassConserve_latlon015.PMID.20180701T20180801.nc'
SurfP_filename = casename+'.cam.h2.MassConserve_latlon015.PS.20180701T20180801.nc'
T_filename = casename+'.cam.h2.MassConserve_latlon015.T.20180701T20180801.nc'

if varname=='HCHO':
    varname_filename = casename+'.cam.h2.MassConserve_latlon015.CH2O.20180701T20180801.nc'
elif varname=='NO2':
    varname_filename = casename+'.cam.h2.MassConserve_latlon015.NO2.20180701T20180801.nc'

#================================================================================================
#### Module import ###
import os
import glob

import pandas as pd
import numpy as np # for array manipulation and basic scientific calculation
import xarray as xr # To read NetCDF files
from scipy.interpolate import griddata # Simple regridding
from netCDF4 import Dataset # To write NetCDF files

# interpolation
from scipy.interpolate import griddata # Simple regridding

import warnings
# Ignore a specific warning
warnings.filterwarnings('ignore', message='Some warning message')
# Ignore the FutureWarning caused by iteritems()
warnings.filterwarnings("ignore", category=FutureWarning, message=".*iteritems.*")

#================================================================================================
### Self-defined functions
# Using the following functions that re-calculate TROPOMI AK after vertically interpolated to MUSICA
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../functions/'))
from func_MUSICA_DefineRegion import *
from func_VerticalRegrid_TROPOMIAK_toMUSICAlevs_015latlon import *

# Pressure functions
def Pres_TO_Height(air_Pres,P_surf):
    """This function calculates altitude (m) based on pressure (Pa), given surface pressure (P_surf)"""
    # Use scale height H ~ 8.5km
    H = 8.5 # km
    # P0 = 101325 #Pa
    heighti = -H*np.log(air_Pres/P_surf)*1e3 # convert to m
    return heighti

#================================================================================================
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
### Loop through all regions
for fileregion in fileregion_ls:
    # Outputfile name and path
    SavePath = MUSICATropVCDSave_diri+casename+'.TropVCD1330LT_withTROPOMIAK.'+varname+'.'+fileregion+'.'+startdate+'T'+enddate+'.nc'
    if os.path.exists(SavePath):
        print(f"The file {SavePath} exists.")
    else:
        print(f'Start to process for {varname}, over {fileregion}')
        # Get the dimensions for a region
        lon_right,lat_bot,lon_left,lat_up = latlonbound(fileregion)

        regionivar_ds = xr.open_dataset( regridByVar_diri+varname_filename ).sel(lat=slice(lat_bot,lat_up),lon=slice(lon_right,lon_left))
        # get the horizontal resolution grid
        regionilats = regionivar_ds.lat
        regionilons = regionivar_ds.lon

        #================================================================================================
        # Size of final dataset to store tropVCD (end product)
        TropVCD_ar = np.zeros([len(TROPOMI_days),len(regionilats),len(regionilons)])
        TropVCD_ar[:,:,:] = np.nan

        #================================================================================================
        # Loop through date, lat, lon
        for dateidx in range(len(TROPOMI_days)):
            datei = TROPOMI_days[dateidx]
            dateifmt2 = TRROPOMIdate_onMUSICAdate[datei]  
            # print('Process %s over %s'%(dateifmt2,fileregion))
            #================================================================================================
            # Get the MUSICA Regrid data and select for the given date, approx at 13:30 PM
            # Approx at TROPOMI overpass time for datei
            averageTimeStart,averageTimeEnd = timezone_approx1330LT(fileregion)
            sel_times = [np.datetime64(dateifmt2+averageTimeStart),np.datetime64(dateifmt2+averageTimeEnd)]
            # Read in: select for the given region & approx two hours at 1 and 2 PM, then take the average
            regioniTROP_P_1330LT_ds = xr.open_dataset( regridByVar_diri+TROP_P_filename ).sel(lat=slice(lat_bot,lat_up),lon=slice(lon_right,lon_left),time=sel_times).mean(dim='time')  
            regioniPMID_1330LT_ds = xr.open_dataset( regridByVar_diri+PMID_filename ).sel(lat=slice(lat_bot,lat_up),lon=slice(lon_right,lon_left),time=sel_times).mean(dim='time')
            regioniSurfP_1330LT_ds = xr.open_dataset( regridByVar_diri+SurfP_filename ).sel(lat=slice(lat_bot,lat_up),lon=slice(lon_right,lon_left),time=sel_times).mean(dim='time')  
            regioniT_1330LT_ds = xr.open_dataset( regridByVar_diri+T_filename ).sel(lat=slice(lat_bot,lat_up),lon=slice(lon_right,lon_left),time=sel_times).mean(dim='time') 
            regionivar_1330LT_ds = regionivar_ds.sel(time=sel_times).mean(dim='time')

            #================================================================================================
            # get the vertically interpolated TROPOMI AK the given varname
            vertically_gridded_AK_da = func_VerticalRegrid_TROPOMIAK_toMUSICAlevs_015latlon(casename, datei, fileregion, varname)

            ## Calculate tropospheric VCD
            # loop through the (lat,lon)
            for latidx in range(len(regionilats)):
                lati = regionilats[latidx]
                for lonidx in range(len(regionilons)):
                    loni = regionilons[lonidx]
                    # print(latidx,lonidx)
                    tropopauseP_Pa = float(regioniTROP_P_1330LT_ds.TROP_P.sel(lat=lati,lon=loni).values) # vertical pressure in Pa
                    SurfP_Pa = float(regioniSurfP_1330LT_ds.PS.sel(lat=lati,lon=loni).values) # find surface pressure in Pa
                    allPMID_Pa = regioniPMID_1330LT_ds.PMID.sel(lat=lati,lon=loni) # vertical pressure in Pa
                    # get the index of lev of PMID that is lower than tropopause (Pressure > Tropopause pressure)
                    # Find the indices where allPMID_Pa is greater than tropopauseP_Pa
                    indices = np.where(allPMID_Pa.values >= tropopauseP_Pa)[0] # to an array

                    ### Get everything below the tropopause (ordered from tropopause to the surface)
                    PMID_latiloni_Pa = allPMID_Pa.isel(lev=slice(indices[0],indices[-1]+1)).values
                    T_latiloni_K = regioniT_1330LT_ds.isel(lev=slice(indices[0],indices[-1]+1)).sel(lat=lati,lon=loni).T.values
                    tropAKs = vertically_gridded_AK_da.isel(MUSICA_layer=slice(indices[0],indices[-1]+1)).sel(lat=lati,lon=loni).values
                    if varname=='HCHO':
                        # adjust different naming approach
                        MUSICAvar_latiloni = regionivar_1330LT_ds.sel(lat=lati,lon=loni).isel(lev=slice(indices[0],indices[-1]+1))['CH2O'].values
                    else:
                        MUSICAvar_latiloni = regionivar_1330LT_ds.sel(lat=lati,lon=loni).isel(lev=slice(indices[0],indices[-1]+1))[varname].values

                    troplayernum = len(MUSICAvar_latiloni) 
                    # # calculate for the given variable, output sum up
                    latloni_AKadjusted_var_molec_cm2 = []

                    #####
                    # # Loop through the vertical layers within the troposphere, 
                    for LAYidx, reversedLAYidx in enumerate(range(troplayernum - 1, -1, -1)):
                        # backwards from surface to Tropopause (currently everything is order from TOA to surface)
                        MUSICA_levi_PMID = PMID_latiloni_Pa[reversedLAYidx] # Pa
                        MUSICA_levi_TEMP = T_latiloni_K[reversedLAYidx] # K
                        AK_levi = tropAKs[reversedLAYidx] # unitless
                        var_moleFraci = MUSICAvar_latiloni[reversedLAYidx] # mole var/mole air

                        if reversedLAYidx==(troplayernum-1):
                        # If for the near-surface level, use the height to the surface
                            MUSICA_levi_H_cm = Pres_TO_Height(MUSICA_levi_PMID,SurfP_Pa)*(1e2)
                        else:
                            # if not at the near-surface level, get the pressure of the layer below as well and calculate the diff
                            MUSICA_levibelow_PMID = PMID_latiloni_Pa[reversedLAYidx+1] # Pa float(MUSICA_levibelow_ds.PMID.values)
                            MUSICA_levi_H_cm = Pres_TO_Height(MUSICA_levi_PMID,SurfP_Pa)*(1e2)-Pres_TO_Height(MUSICA_levibelow_PMID,SurfP_Pa)*(1e2)

                        # Calculate the sum of all cells across a layer (unit originally in mol/mol)
                        Avog = 6.0221409e23
                        m3Tocm3 = 1e6
                        R = 8.314462 
                        DryAirMass = 0.029 # kg/mol
                        rho_airi = (MUSICA_levi_PMID*DryAirMass)/(R*MUSICA_levi_TEMP)

                        levi_var_molec_cm2 = var_moleFraci*Avog*(1/DryAirMass)*rho_airi*(1/m3Tocm3)*MUSICA_levi_H_cm

                        AKadjusted_levi_var_molec_cm2 = levi_var_molec_cm2*AK_levi
                        # add to the array to calculate the sum of all levels
                        latloni_AKadjusted_var_molec_cm2.append(AKadjusted_levi_var_molec_cm2)

                    # (lat,lon) add tropVCD summed up for all layers, put to the corresponding lat, lon grid
                    TropVCD_ar[dateidx,latidx,lonidx] = np.nansum(np.array(latloni_AKadjusted_var_molec_cm2))

        #================================================================================================
        ## Write to Xr Dataset and save
        import datetime
        today_date = datetime.date.today()

        # write this variable to a dataarray
        MUSICA_TropVCD_ds = xr.Dataset(
                            {
                                varname+'_TropVCD': (['time','lat', 'lon'], TropVCD_ar),
                            },
                            coords={'time': ('time', MUSICA_days), 
                                    'lat': ('lat', regionilats.values), 
                                    'lon': ('lon', regionilons.values)},
                            attrs={'process date': str(today_date.strftime("%Y-%m-%d")),
                                   'description': 'MUSICA VCD adjusted using Averaging Kernels from TROPOMI',
                                  'region': fileregion,},
                            )

        # Add the attributes of coordinates
        for var_name in ['lat','lon']:    
            # Copy attributes from source to target variable
            for attr_name, attr_value in regionivar_ds[var_name].attrs.items():
                MUSICA_TropVCD_ds[var_name].attrs[attr_name] = attr_value

        MUSICA_TropVCD_ds[varname+'_TropVCD'].attrs['longname'] = 'tropospheric vertical column density of '+varname
        MUSICA_TropVCD_ds[varname+'_TropVCD'].attrs['unit'] = 'molec/cm2'

        #================================================================================================
        # save to .nc
        MUSICA_TropVCD_ds.to_netcdf(SavePath)
        print('Save to:',SavePath)
