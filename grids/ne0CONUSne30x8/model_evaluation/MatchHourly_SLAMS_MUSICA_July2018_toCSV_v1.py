'''
This script is used to get matched MUSICA outputs in ne0CONUSne30x8 horizontal grid with SLAMS/AQS hourly measurements

MODIFICATION HISTORY:
    Madankui Tao, 18, DEC, 2023: VERSION 1.0
    - Initial version
'''
#================================================================================================
# Specified variables
fileregion_ls = ['Northeast']
# fileregion_ls = ['WestCoast','Mountain','Midwest','Southwest','Southeast','Northeast']
#Options: 'WestCoast','Mountain','Midwest','Southwest','Southeast','Northeast'

### Date
startdate = "2018-07-01" #"2018-06-30"
enddate = "2018-08-01" #"2018-08-03"
# filetime = f'{startdate}T{enddate}'

# casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1base01'
# casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1basehourlyNEI2017'
casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1baseJulyMeanNEI2017'

# ============================================================
# USER CONFIGURATION — set these paths before running
# ============================================================
AQS2018_diri = ''           # Path to AQS data directory for year 2018 (contains hourly zip files)
AQS_diri = ''               # Path to base AQS data directory (contains parameters.csv)
Matched_diri = ''           # Path to directory for output matched CSV files
MUSICA_index_diri = ''      # Path to directory containing MUSICA column-index CSV files
svante_archive = ''         # Path to your CESM archive directory
# ============================================================

# Where to store the matched files | since calculation takes too long, store for each date first
# Matched_diri = f'.../{casename}/hourly/ByEachDay/'

#------------------------------------------------------------------------------
### For AQS
hourly_varlist = ['O3','CO','NO2',]
# hourly_varlist = ['CO','SO2','NO2','O3','PM25']
# 'LC25' for PM2.5 - Local Conditions (88101)

# for hourly measurement
valcolname = 'Sample Measurement'

# Name used by AQS
AQSvarindex_dic = {'CO':'CO','SO2':'SO2',
                   'PM25':'LC25',
                   'O3':'O3',
                   'NO2':'NO2',
                   'HCHO':'Formaldehyde','ISOP':'Isoprene','Benzene':'Benzene'}

#------------------------------------------------------------------------------
### For MUSICA
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

# "Formaldehyde" in parameters_df["Parameter"].unique()
# HCHO_code = ParaNameCode_dic["Formaldehyde"]
# ISOP_code = ParaNameCode_dic["Isoprene"]    

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
for varname in hourly_varlist:
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
            AQS_filename = 'hourly_'+str(ParaAbrrCode_dic[varname])+'_2018.zip'
            AQS_zf = ZipFile(AQS2018_diri+AQS_filename)
            # see what's in the zip file
            Inside_filenames = deflatten(AQS_zf.namelist())
            # read in the csv file in the zip file ,
            fields = ['State Code','County Code','Site Num','Parameter Code','POC','Latitude','Longitude','Date Local','Time Local',valcolname,'Units of Measure']
            AQS_df = pd.read_csv(AQS_zf.open(Inside_filenames[0]),skipinitialspace=True, index_col=False, usecols=fields)
            # # Add a new column combining State Code-County Code-Site Num
            AQS_df = AQS_df.astype({col: str for col in ["State Code", "County Code", "Site Num", "Parameter Code"]})
            AQS_df["MonitorID"] = AQS_df[["State Code", "County Code", "Site Num","Parameter Code"]].agg('-'.join, axis=1) 
            # Simplify the df
            AQS_df = AQS_df[['MonitorID','Latitude','Longitude','Date Local','Time Local',valcolname,'Units of Measure']]
            # Select only for between start and end date
            AQS_selTime_df = AQS_df[(AQS_df['Date Local']>=startdate)&(AQS_df['Date Local']<=enddate)]
            # Select for the given region
            regioni_AQS_selTime_df = AQS_selTime_df[(AQS_selTime_df['Latitude']>=lat_bot)&(AQS_selTime_df['Latitude']<=lat_up)
                                             & (AQS_selTime_df['Longitude']>=lon_right)&(AQS_selTime_df['Longitude']<=lon_left)]
            
            #================================================================================================
            ### hourly (h2) MUSICA added as a new column
            # svante_archive is set in the USER CONFIGURATION block at the top of this file
            RunPath = svante_archive+casename+'/atm/hist/'
            
            # Loop through each MonitorID, using it as the key to get the MUSICA column index
            # Get the matched MUSICA outputs and to the same unit
            import pytz
            from datetime import datetime

            # Match with surface simulations save to an array (use array instead of list for efficiency)
            MUSICA_hourlyVar_ar = np.zeros(regioni_AQS_selTime_df.shape[0])
            MUSICA_hourlyVar_ar[:] = np.nan

            # Get the correct unit conversion
            if regioni_AQS_selTime_df['Units of Measure'].values[0]=='Parts per billion':
                scalefactor=1e9
            elif regioni_AQS_selTime_df['Units of Measure'].values[0]=='Parts per million':
                scalefactor=1e6
            elif regioni_AQS_selTime_df['Units of Measure'].values[0]=='Micrograms/cubic meter (LC)':
                # For PM2.5, convert from kg/m3 to Micrograms/m3
                scalefactor=1e9
            else:
                print('Stop and check the unit')

            print('Processing number of rows: ',str(regioni_AQS_selTime_df.shape[0]))
            for idx in range(regioni_AQS_selTime_df.shape[0]):
                MonitorIDi = regioni_AQS_selTime_df.MonitorID.values[idx]
                LT_Datei = regioni_AQS_selTime_df['Date Local'].values[idx]
                LT_Timei = regioni_AQS_selTime_df['Time Local'].values[idx]

                # Create a datetime object from the local time and set the timezone
                local_datetime = datetime.strptime(LT_Datei + ' ' + LT_Timei, '%Y-%m-%d %H:%M')
                # Calculate the final UTC time by adding the offset
                final_utc_time = local_datetime + timedelta(hours=timezone_offset)
                # Convert the final UTC time to numpy datetime64 format
                final_utc_time_dt64 = np.datetime64(final_utc_time)

                ### Get MUSICA value
                MUSICA_filedatei = final_utc_time.strftime('%Y-%m-%d')
                # Read in the hourly mean file
                datei_MUSICAfilePath = f'{RunPath}{casename}.cam.h2.{MUSICA_filedatei}-03600.nc'
                if os.path.exists(datei_MUSICAfilePath)==False:
                    # print(f"MUSICA file on {MUSICA_filedatei} missing.")
                    # MUSICA_hourlyVar.append(np.nan) # keep it nan
                    pass
                else:
                    # Proceed to extract MUSICA
                    datei_MUSICA_ds = xr.open_dataset(datei_MUSICAfilePath)
                    # if both time and MonitorIDi are valid keys to get MUSICA val
                    if (final_utc_time_dt64 in datei_MUSICA_ds.time.values) & (MonitorIDi in MUSICAcolIndex_dic):
                        # avoid the issue of missing specific time
                        # select for the given time, location at the surface
                        selsurf_ds = datei_MUSICA_ds.sel(time=final_utc_time_dt64,ncol=int(MUSICAcolIndex_dic[MonitorIDi])).isel(lev=-1)
                        vari_MUSICA_val = float(selsurf_ds[MUSICA_varname_dic[varname]].values)*scalefactor
                        # add to the array
                        MUSICA_hourlyVar_ar[idx] = vari_MUSICA_val
                        # MUSICA_hourlyVar.append(vari_MUSICA_val)

                    else:
                        # print(f"Missing time {final_utc_time_dt64}")
                        # MUSICA_hourlyVar.append(np.nan) # keep it nan
                        pass

            # Add the array as a column in the df
            regioni_AQS_selTime_df['MUSICAv0_approxLThourly'] = MUSICA_hourlyVar_ar
            # drop the nan values
            nona_regioni_AQS_selTime_df = regioni_AQS_selTime_df.dropna()
            print('Matched df size',nona_regioni_AQS_selTime_df.shape)            
            #================================================================================================
            # Save to .csv file
            nona_regioni_AQS_selTime_df.to_csv(Savefile_path, index=False)
            print(f'Saved to: {Savefile_path}')

#================================================================================================
