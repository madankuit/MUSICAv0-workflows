'''
This script is used to get matched Daily mean MUSICA outputs (converted to LT) in ne0CONUSne30x8 horizontal grid with SLAMS/AQS hourly measurements
*MUSICA saved the previous mean to the current one

MODIFICATION HISTORY:
    M. Tao, 19, DEC, 2023: VERSION 1.0
    - Initial version
'''
#================================================================================================
# Specified variables
# fileregion_ls = ['Northeast']
fileregion_ls = ['WestCoast','Mountain','Midwest','Southwest','Southeast','Northeast']

### Date
startdate = "2018-07-01" #"2018-06-30"
enddate = "2018-08-01" #"2018-08-03"
filetime = f'{startdate}T{enddate}'

# casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1base01'
# casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1basehourlyNEI2017'
casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1basehourlyNOotherJulyMeanNEI2017'
# casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1baseJulyMeanNEI2017'
### with 6-hr nudging
# casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.6HrNudgeTS1hourlyNEI2017'

# ============================================================
# USER CONFIGURATION — set these paths before running
# ============================================================
AQS2018_diri = ''           # Path to AQS data directory for year 2018 (contains daily zip files)
AQS_diri = ''               # Path to base AQS data directory (contains parameters.csv)
Matched_diri = ''           # Path to directory for output matched CSV files
MUSICA_index_diri = ''      # Path to directory containing MUSICA column-index CSV files
svante_archive = ''         # Path to your CESM archive directory
# ============================================================

# Where to store the matched files
# Matched_diri = f'.../{casename}/daily/'

#------------------------------------------------------------------------------
### For AQS
# daily_varlist = ['NO2']
daily_varlist = ['NO2','O3','CO','SO2','PM25']
# 'LC25' for PM2.5 - Local Conditions (88101)

# for daily measurement
valcolname = 'Arithmetic Mean'

# Name used by AQS
AQSvarindex_dic = {'CO':'CO','SO2':'SO2',
                   'PM25':'LC25',
                   'O3':'O3',
                   'NO2':'NO2',
                   'HCHO':'Formaldehyde','ISOP':'Isoprene','Benzene':'Benzene'}

#------------------------------------------------------------------------------
### For MUSICA
RunPath = svante_archive+casename+'/atm/hist/'
# MUSICA_index_diri is set in the USER CONFIGURATION block above

# Name used in MUSICA outputs
MUSICA_varname_dic = {'CO':'CO','SO2':'SO2',
                   'PM25':'PM25',
                   'O3':'O3',
                   'NO2':'NO2',
                   'HCHO':'CH2O','ISOP':'ISOP','Benzene':'BENZENE'}

#================================================================================================
# ### functions import ###

# basics
import os
import os.path
from os.path import exists
from os import path
import sys  
from io import BytesIO
import requests
from zipfile import ZipFile
import json

# numerical
import math
import random
import numpy as np
import numpy.ma as ma
import pandas as pd
import netCDF4 as nc4
from netCDF4 import Dataset
import xarray as xr

# my functions
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../functions/'))
from func_MUSICA_DefineRegion import *

# Grid
# Read SCRIP file that has grid information needed to plot values on a map
SCRIP_CONUS = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../grid_files/ne0CONUS_ne30x8_np4_SCRIP.nc')
SCRIP_ne30 = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../grid_files/ne30np4_091226_pentagons.nc')

#================================================================================================
### """Functions for AQS"""
from zipfile import ZipFile

def deflatten(names):
    """This function create a sort of structure from the file list of a zip file"""
    names.sort(key=lambda name:len(name.split('/')))
    deflattened = []
    while len(names) > 0:
        name = names[0]
        if name[-1] == '/':
            subnames = [subname[len(name):] for subname in names if subname.startswith(name) and subname != name]
            for subname in subnames:
                names.remove(name+subname)
            deflattened.append((name, deflatten(subnames)))
        else:
            deflattened.append(name)
        names.remove(name)
    return deflattened

# get parameter code
parameters_df = pd.read_csv(AQS_diri+"parameters.csv")
ParaAbrrCode_dic = dict(zip(parameters_df["Parameter Abbreviation"].values, parameters_df["Parameter Code"].values))
ParaAbrrName_dic = dict(zip(parameters_df["Parameter Abbreviation"].values, parameters_df["Parameter"].values))
ParaAbrrUnit_dic = dict(zip(parameters_df["Parameter Abbreviation"].values, parameters_df["Standard Units"].values))

#================================================================================================
# Date format dictionary for processing dates
# Find the list of dates
from datetime import datetime, timedelta
from pytz import timezone

fmt1_days = []
fmt2_days = []
import datetime
Read_start = datetime.datetime.strptime(startdate, "%Y-%m-%d")
Read_end_plus1 = datetime.datetime.strptime(enddate, "%Y-%m-%d")
date_generated = [Read_start + timedelta(days=x) for x in range(0, (Read_end_plus1-Read_start).days)]

# In the format of filename date
for date in date_generated:
    fmt1_days.append(date.strftime("%Y%m%d"))
    fmt2_days.append(date.strftime("%Y-%m-%d"))
    
# in a dictionary to adjust for different date format
date_fmt_convertdic = dict(zip(fmt1_days,fmt2_days))

#================================================================================================
#================================================================================================
### Loop through variables
for varname in daily_varlist:
    # print(varname)
    #-----------------------------------------------------
    # Get the macthed MUSICA column index
    vari_colidxfile_path = MUSICA_index_diri+'MUSICA0_hourlyAQS_'+AQSvarindex_dic[varname]+'_Monitor_colidx.csv'
    allvari_colidx_df = pd.read_csv(vari_colidxfile_path)
    # Select only for sites that have good matching with the model
    vari_colidx_df = allvari_colidx_df[allvari_colidx_df['MUSICA0_colIndex']!='Find None']

    # Use MonitorID as the key for MUSICA0_colIndex
    MUSICAcolIndex_dic = dict(zip(vari_colidx_df["MonitorID"].values, vari_colidx_df["MUSICA0_colIndex"].values))
    # print(varname,'monitor number in the CONUS:',vari_colidx_df.shape[0])
    
    ### Loop through all regions
    for fileregion in fileregion_ls:
        # Outputfile name and path
        Savefile_path = f'{Matched_diri}{casename}.MatchedSLAMS.{fileregion}.{varname}.{filetime}.csv'
        if os.path.exists(Savefile_path):
            print(f"{Savefile_path} already exists.")
        else:
            print(f"Start processing for {casename}, {varname} over  {fileregion} | {filetime}")
            # proceed 
            lon_right,lat_bot,lon_left,lat_up = latlonbound(fileregion)
            # Get MUSICA time (UTC)
            timezone_offset = UTCtimezone_offset_summer(fileregion)  # UTC-X
            #================================================================================================
            ### AQS as the base dataframe
            AQS_filename = 'daily_'+str(ParaAbrrCode_dic[AQSvarindex_dic[varname]])+'_2018.zip'
            AQS_zf = ZipFile(AQS2018_diri+AQS_filename)
            # see what's in the zip file
            Inside_filenames = deflatten(AQS_zf.namelist())
            # read in the csv file in the zip file ,
            fields = ['State Code','County Code','Site Num','Parameter Code','POC','Latitude','Longitude','Date Local',valcolname,'Units of Measure']
            AQS_df = pd.read_csv(AQS_zf.open(Inside_filenames[0]),skipinitialspace=True, index_col=False, usecols=fields)
            # # Add a new column combining State Code-County Code-Site Num
            AQS_df = AQS_df.astype({col: str for col in ["State Code", "County Code", "Site Num", "Parameter Code"]})
            AQS_df["MonitorID"] = AQS_df[["State Code", "County Code", "Site Num","Parameter Code"]].agg('-'.join, axis=1) 
            # Simplify the df
            AQS_df = AQS_df[['MonitorID','Latitude','Longitude','Date Local',valcolname,'Units of Measure']]
            # Select only for between start and end date
            AQS_selTime_df = AQS_df[(AQS_df['Date Local']>=startdate)&(AQS_df['Date Local']<=enddate)]
            # Select for the given region
            regioni_AQS_selTime_df = AQS_selTime_df[(AQS_selTime_df['Latitude']>=lat_bot)&(AQS_selTime_df['Latitude']<=lat_up)
                                                     & (AQS_selTime_df['Longitude']>=lon_right)&(AQS_selTime_df['Longitude']<=lon_left)]
            
            ndfrows = regioni_AQS_selTime_df.shape[0]
            print('Processing number of rows: ',str(ndfrows))
            
            #================================================================================================
            # Get the correct unit conversion
            
            AQSunit = regioni_AQS_selTime_df['Units of Measure'].values[0]
            
            if AQSunit=='Parts per billion':
                scalefactor=1e9
            elif AQSunit=='Parts per million':
                scalefactor=1e6
            elif AQSunit=='Micrograms/cubic meter (LC)':
                # For PM2.5, convert from kg/m3 (MUSICA) to Micrograms/m3
                scalefactor=1e9
            else:
                print('Stop and check the unit')
            # print(varname,AQSunit,scalefactor)
                
            #================================================================================================
            ### Local Time Daily Mean using hourly (h2) MUSICA; added as a new column
            
            # Loop through each MonitorID, using it as the key to get the MUSICA column index
            # Get the matched MUSICA outputs and to the same unit
            import pytz
            from datetime import datetime

            # Match with surface simulations save to an array (use array instead of list for efficiency)
            MUSICA_dailyVar_ar = np.zeros([ndfrows])
            MUSICA_dailyVar_ar[:] = np.nan

            for idx in range(ndfrows):
                MonitorIDi = regioni_AQS_selTime_df.MonitorID.values[idx]
                LT_Datei = regioni_AQS_selTime_df['Date Local'].values[idx]

                # find the start and end times in LT from this day (T00:00 to T:23:00)
                Datei00_LTdatetime = datetime.strptime(LT_Datei + ' 00:00', '%Y-%m-%d %H:%M')
                Datei23_LTdatetime = datetime.strptime(LT_Datei + ' 23:00', '%Y-%m-%d %H:%M')
                # Calculate the final UTC time by adding the offset
                Datei00_utc_time = Datei00_LTdatetime + timedelta(hours=timezone_offset)
                Datei23_utc_time = Datei23_LTdatetime + timedelta(hours=timezone_offset)
                # Convert the final UTC time to numpy datetime64 format
                Datei00_utc_time_dt64 = np.datetime64(Datei00_utc_time)
                Datei23_utc_time_dt64 = np.datetime64(Datei23_utc_time)

                # Read in both date h2 files
                Datei00_MUSICA_filedatei = Datei00_utc_time.strftime('%Y-%m-%d')
                Datei23_MUSICA_filedatei = Datei23_utc_time.strftime('%Y-%m-%d')

                # Check if MonitorIDi is a valid key to get MUSICA val
                if MonitorIDi in MUSICAcolIndex_dic:
                    ### Get MUSICA value if the column index is valid
                    # verify files exist
                    Datei00_MUSICAfilePath = f'{RunPath}{casename}.cam.h2.{Datei00_MUSICA_filedatei}-03600.nc'
                    Datei23_MUSICAfilePath = f'{RunPath}{casename}.cam.h2.{Datei23_MUSICA_filedatei}-03600.nc'
                    if (os.path.exists(Datei00_MUSICAfilePath)==True) and (os.path.exists(Datei23_MUSICAfilePath)==True):
                        # Proceed to extract MUSICA variable at the surface for the given location
                        ## For the ones falling on the first UTC time day
                        # Logic: surface (lev=-1); time; location using ncol; 
                        Datei00_MUSICA_da = xr.open_dataset(Datei00_MUSICAfilePath).isel(lev=-1).sel(time=slice(Datei00_utc_time_dt64,None),ncol=int(MUSICAcolIndex_dic[MonitorIDi]))
                        # and given variable
                        Datei00_var_da = Datei00_MUSICA_da[MUSICA_varname_dic[varname]]#.values
                        ## For the ones falling on the second UTC time day
                        Datei23_MUSICA_da = xr.open_dataset(Datei23_MUSICAfilePath).isel(lev=-1).sel(time=slice(None,Datei23_utc_time_dt64),ncol=int(MUSICAcolIndex_dic[MonitorIDi]))
                        Datei23_var_da = Datei23_MUSICA_da[MUSICA_varname_dic[varname]]
                        # Merge
                        Datei_var_da = xr.concat([Datei00_var_da, Datei23_var_da], dim="time")
                        vari_MUSICA_val = np.nanmean(Datei_var_da.values)*scalefactor
                        # add to the array
                        MUSICA_dailyVar_ar[idx] = vari_MUSICA_val  
                    else:
                        print(f"MUSICA file on {Datei00_MUSICAfilePath} or {Datei23_MUSICAfilePath} missing.")
                        # MUSICA_hourlyVar.append(np.nan) # keep it nan

            # print('MUSICA array size',len(MUSICA_dailyVar_ar))
            # Add the array as a column in the df
            regioni_AQS_selTime_df['MUSICAv0_approxLTdaily'] = MUSICA_dailyVar_ar
            # drop the nan values
            nona_regioni_AQS_selTime_df = regioni_AQS_selTime_df.dropna()
            print('Matched df size',nona_regioni_AQS_selTime_df.shape)

            #================================================================================================
            # Save to .csv file
            nona_regioni_AQS_selTime_df.to_csv(Savefile_path, index=False)
            print(f'Saved to: {Savefile_path}')

# #================================================================================================
