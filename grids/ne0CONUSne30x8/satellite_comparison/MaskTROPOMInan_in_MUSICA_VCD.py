# Go through each date and masked out TROPOMI nan values in MUSICA-tropVCD

# Input: Grided L3 HiRv2 TROPOMI & MUSICA-tropVCD 
# Output: Masked MUSICA-tropVCD where TROPOMI has nan

#=====Processing period (can be changed)======

startDate = "2018-07-01"
endDate = "2018-08-01" 

varname = 'CO' # options: HCHO, NO2 (use tropospheric VCD); CO (use total VCD)

filetime = f'{startDate}T{endDate}'

fileregion_ls = ['WestCoast','Mountain','Midwest','Southwest','Southeast','Northeast']

casename_ls = [
                # 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1base01',
                # 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1basehourlyNEI2017',
                # 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1baseJulyMeanNEI2017',
    'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1basehourlyNOotherJulyMeanNEI2017',
    # 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.6HrNudgeTS1hourlyNEI2017',
              ]

# TROPOMI data info
fileheader_dic = {'HCHO':'S5P_RPRO_L2__HCHO___',
               'NO2':'S5P_RPRO_L2__NO2____',
               'CO':'S5P_RPRO_L2__CO_____',
              }

L3ProductFileOut_diri_dic = {'HCHO':'/net/fs09/d0/taoma528/Datasets/S5P_L2__HCHO___HiR_2/regrid_2018_MUSICA015/',
                               'NO2':'/net/fs09/d0/taoma528/Datasets/S5P_L2__NO2____HiR_2/regrid_2018_MUSICA015/',
                               'CO':'/net/fs09/d0/taoma528/Datasets/S5P_L2__CO_____HiR_2/regrid_2018_MUSICA015/',
                              }

#--------------------------------------------------------------------------
#=====Dependent library & functions======
#--------------------------------------------------------------------------
# basics
import sys  
from io import BytesIO
import requests
from zipfile import ZipFile
import os.path
from os import path

# numerical
import math
import numpy as np
import numpy.ma as ma
import pandas as pd
import h5py
import netCDF4 # as nc4
from netCDF4 import Dataset
import xarray as xr
xr.set_options(display_style="html")
from sklearn.datasets import make_regression
from sklearn.linear_model import LinearRegression

# datetime
import datetime
from datetime import datetime, timedelta
from pytz import timezone

# set to ignore warnings
import warnings
warnings.filterwarnings("ignore")

import sys
sys.path.insert(0,'/home/taoma528/Scripts/CESM_analysis/functions/')
from func_MUSICA_DefineRegion import *

#=====Get the processing dates======
#--------------------------------------------------------------------------
GivenDays = []
import datetime
# from datetime import datetime
start = datetime.datetime.strptime(startDate, "%Y-%m-%d")
end = datetime.datetime.strptime(endDate, "%Y-%m-%d")
date_generated = [start + datetime.timedelta(days=x) for x in range(0, (end-start).days)] # exclude end date
# next_date = end + datetime.timedelta(days=1)
# date_generated = [start + datetime.timedelta(days=x) for x in range(0, (next_date-start).days)] # include end date

for date in date_generated:
    GivenDays.append(date.strftime("%Y-%m-%d"))
print('Read in %i Day(s)'%len(GivenDays))

#================================================================================================
#================================================================================================
if varname in ['NO2','HCHO']:
    #================================================================================================
    #### Tropospheric VCD
    ### Loop through all cases and regions
    for casename in casename_ls:
        # MUSICA-TropVCD directory 
        MUSICATropVCD_diri = f'/net/fs09/d0/taoma528/CESM22/Regridded_MUSICA_Output/2018_1330LT_TROPOMIcomp/MassConserve_latlon015_MUSICAoutput/{casename}/TropVCD_approx1330LT/'

        for fileregion in fileregion_ls:
            lon_right,lat_bot,lon_left,lat_up = latlonbound(fileregion)
            # Outputfile name and path
            SavePath = f'{MUSICATropVCD_diri}MaskedTROPOMInan_{casename}.TropVCD1330LT_withTROPOMIAK.{varname}.{fileregion}.{filetime}.nc'

            if os.path.exists(SavePath):
                print(f"The file {SavePath} exists.")
            else:
                print(f'Start to process for {varname}, over {fileregion}')
                #--------------------------------------------------------------------------
                # Read in MUSICA for this file region
                regioni_MUSICA_filename = f'{MUSICATropVCD_diri}{casename}.TropVCD1330LT_withTROPOMIAK.{varname}.{fileregion}.{filetime}.nc'
                regioni_MUSICA_ds = xr.open_dataset(regioni_MUSICA_filename)

                ## Loop through to process each day
                #--------------------------------------------------------------------------
                for DayIdx in range(len(GivenDays)):
                    datei = GivenDays[DayIdx]
                    TROPOMIfmt_datei = ('').join(datei.split('-'))
                    # print('Start process for :',GivenDay)
                    #--------------------------------------------------------------------------
                    # Find matching TROPOMI data
                    TROPOMI_tropVCD_filename = L3ProductFileOut_diri_dic[varname]+fileheader_dic[varname]+'Global_Regrid015deg_'+TROPOMIfmt_datei+'.nc'
                    datei_regioni_TROPOMI_tropVCD_ds = xr.open_dataset(TROPOMI_tropVCD_filename).sel(lat=slice(lat_bot,lat_up),lon=slice(lon_right,lon_left),time=datei)
                    #--------------------------------------------------------------------------
                    # Get MUSICA on this date
                    datei_regioni_MUSICA_ds = regioni_MUSICA_ds.sel(time=datei,lat=slice(lat_bot,lat_up),lon=slice(lon_right,lon_left)) #,lat=slice(lat_bot,lat_up),lon=slice(lon_right,lon_left)

                    # Mask out where TROPOMI has NaN
                    datei_regioni_MUSICA_ds[varname+'_TropVCD'] = datei_regioni_MUSICA_ds[varname+'_TropVCD'].where(~np.isnan(datei_regioni_TROPOMI_tropVCD_ds['TROPOMI_'+varname+'_TropVCD']), np.nan)
                    # datei_regioni_MUSICA_ds[varname+'_TropVCDPrecision'] = datei_regioni_MUSICA_ds[varname+'_TropVCDPrecision'].where(~np.isnan(datei_regioni_TROPOMI_tropVCD_ds['TROPOMI_'+varname+'_TropVCD']), np.nan)

                    # if selecting var
                    selected_ds = datei_regioni_MUSICA_ds[varname + '_TropVCD'] #datei_regioni_MUSICA_ds.sel(variable=[varname+'_TropVCD'])

                    if DayIdx==0:
                        merged_ds = selected_ds
                    else:
                        # Merge datasets along the time dimension
                        merged_ds = xr.concat([merged_ds, selected_ds], dim="time")

                # Save the merged_ds
                merged_ds.to_netcdf(SavePath)
                print('Save to:',SavePath)

#================================================================================================
elif varname in ['CO']:
    #================================================================================================
    #### Total VCD
    ### Loop through all cases and regions
    for casename in casename_ls:
        # MUSICA-TotalVCD directory 
        MUSICAVCD_diri = f'/net/fs09/d0/taoma528/CESM22/Regridded_MUSICA_Output/2018_1330LT_TROPOMIcomp/MassConserve_latlon015_MUSICAoutput/{casename}/TotalVCD_approx1330LT/'

        for fileregion in fileregion_ls:
            lon_right,lat_bot,lon_left,lat_up = latlonbound(fileregion)
            # Outputfile name and path
            SavePath = f'{MUSICAVCD_diri}MaskedTROPOMInan_{casename}.TotalVCD1330LT_withTROPOMIAK.{varname}.{fileregion}.{filetime}.nc'

            if os.path.exists(SavePath):
                print(f"The file {SavePath} exists.")
            else:
                print(f'Start to process for {varname}, over {fileregion}')
                #--------------------------------------------------------------------------
                # Read in MUSICA for this file region
                regioni_MUSICA_filename = f'{MUSICAVCD_diri}{casename}.TotalVCD1330LT_withTROPOMIAK.{varname}.{fileregion}.{filetime}.nc'
                regioni_MUSICA_ds = xr.open_dataset(regioni_MUSICA_filename)

                ## Loop through to process each day
                #--------------------------------------------------------------------------
                for DayIdx in range(len(GivenDays)):
                    datei = GivenDays[DayIdx]
                    TROPOMIfmt_datei = ('').join(datei.split('-'))
                    # print('Start process for :',GivenDay)
                    #--------------------------------------------------------------------------
                    # Find matching TROPOMI data
                    TROPOMI_VCD_filename = L3ProductFileOut_diri_dic[varname]+fileheader_dic[varname]+'Global_Regrid015deg_'+TROPOMIfmt_datei+'.nc'
                    datei_regioni_TROPOMI_VCD_ds = xr.open_dataset(TROPOMI_VCD_filename).sel(lat=slice(lat_bot,lat_up),lon=slice(lon_right,lon_left),time=datei)
                    #--------------------------------------------------------------------------
                    # Get MUSICA on this date
                    datei_regioni_MUSICA_ds = regioni_MUSICA_ds.sel(time=datei,lat=slice(lat_bot,lat_up),lon=slice(lon_right,lon_left))

                    # Mask out where TROPOMI has NaN
                    datei_regioni_MUSICA_ds[varname+'_TotalVCD'] = datei_regioni_MUSICA_ds[varname+'_TotalVCD'].where(~np.isnan(datei_regioni_TROPOMI_VCD_ds['TROPOMI_'+varname+'_TotalVCD']), np.nan)
                    # datei_regioni_MUSICA_ds[varname+'_TotalVCDPrecision'] = datei_regioni_MUSICA_ds[varname+'_TotalVCDPrecision'].where(~np.isnan(datei_regioni_TROPOMI_VCD_ds['TROPOMI_'+varname+'_TotalVCD']), np.nan)

                    # if selecting var
                    selected_ds = datei_regioni_MUSICA_ds[varname + '_TotalVCD'] #datei_regioni_MUSICA_ds.sel(variable=[varname+'_TropVCD'])

                    if DayIdx==0:
                        merged_ds = selected_ds
                    else:
                        # Merge datasets along the time dimension
                        merged_ds = xr.concat([merged_ds, selected_ds], dim="time")

                # Save the merged_ds
                merged_ds.to_netcdf(SavePath)
                print('Save to:',SavePath)

