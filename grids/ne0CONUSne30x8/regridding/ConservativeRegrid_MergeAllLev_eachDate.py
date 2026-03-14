'''
This code is designed to merge all vertical layered variables for each day for re-gridded SE outputs of MUSICA-V0 in regular 0.15x0.15 lat and lon grids using mass conservative method

The MUSICAv0 output has already been mass-conservatively regridded into a resolution of 0.15x0.15 lat-lon but in separate files for date and each vertical layer

MODIFICATION HISTORY:
    M. Tao, 4, December, 2023: VERSION 1.0
    - Initial version
    M. Tao, 4, March, 2024: VERSION 1.1
    - Add checks for Merged files
'''

# ============================================================
# USER CONFIGURATION — set these paths before running
# ============================================================
# casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1base01'
# casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1basehourlyNEI2017'
casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1basehourlyNOotherJulyMeanNEI2017'
# casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1baseJulyMeanNEI2017'
### with 6-hr nudging
# casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.6HrNudgeTS1hourlyNEI2017'

# Directory containing per-level regridded h2 files
# e.g. '/path/to/Regridded_MUSICA_Output/.../MassConserve_latlon015_MUSICAoutput/<casename>/h2/'
h2files_diri = ''

# Directory for single-day merged vertical-variable files
# e.g. '/path/to/Regridded_MUSICA_Output/.../MassConserve_latlon015_MUSICAoutput/<casename>/h2_VertiVarsByVarDate/'
singleDay_diri = ''

# Directory for per-variable output files
# e.g. '/path/to/Regridded_MUSICA_Output/.../MassConserve_latlon015_MUSICAoutput/<casename>/h2_ByVar/'
FileOutByVar_dir = ''

# Directory for merged-vertical h2 files (often same as singleDay_diri)
# e.g. '/path/to/Regridded_MUSICA_Output/.../MassConserve_latlon015_MUSICAoutput/<casename>/h2_VertiVarsByVarDate/'
MergedVerticalh2files_diri = ''
# ============================================================

GroupedName = 'vertivars' # 'vertivars' ; '.CO.'
# single species files use .XX. to make sure the file selection was effective

startDate = "2018-07-01"
endDate = "2018-08-01"
endDateExclude = "2018-08-02" # the created datearray exclude this date

### version for case: f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1basehourlyNEI2017
var_dic = {
           'vertivars':['CO'] #,'CH2O','NO2','T','PMID'], #['CO','CH2O','NO2','T','PMID']
          }

### version for case: f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1base01
# var_dic = {'NoVertvars':['TROP_P','PS'],
#            'vertivars':['CH2O','NO2','T','PMID'],
#            '.CO.':['CO'],
#           }

#================================================================================================
### Module import ###
import os
import glob

import pandas as pd
import numpy as np # for array manipulation and basic scientific calculation
import xarray as xr # To read NetCDF files
from netCDF4 import Dataset # To write NetCDF files

# interpolation
from scipy.interpolate import griddata # Simple regridding

import warnings
# Ignore a specific warning
warnings.filterwarnings('ignore', message='Some warning message')
# Ignore the FutureWarning caused by iteritems()
warnings.filterwarnings("ignore", category=FutureWarning, message=".*iteritems.*")

from datetime import datetime
import calendar

import re
def extract_numeric_part(filename):
    # Use regular expression to find all sequences of digits in the filename
    numeric_parts = re.findall(r'\d+', filename)
    
    # Convert the numeric parts to integers and return as a tuple
    numeric_values = [int(part) for part in numeric_parts]
    return tuple(numeric_values)

#=====Get the processing dates======
#--------------------------------------------------------------------------
GivenDays = []
GivenDays_fmt2 = []
import datetime
# from datetime import datetime
start = datetime.datetime.strptime(startDate, "%Y-%m-%d")
end = datetime.datetime.strptime(endDateExclude, "%Y-%m-%d")
date_generated = [start + datetime.timedelta(days=x) for x in range(0, (end-start).days)]

for date in date_generated:
    GivenDays.append(date.strftime("%Y%m%d"))
    GivenDays_fmt2.append(date.strftime("%Y-%m-%d")) 
# print('Read in %i Day(s)'%len(GivenDays))
print(GivenDays[0],GivenDays[-1])

#=====Check if missing files======
#--------------------------------------------------------------------------
# Check and get file names
alldates_files = []
# Read in the grouped names for all times given
for datei in GivenDays_fmt2:
    # Get a list of file paths in the specified directories with the given date
    datei_filelist = glob.glob(os.path.join(h2files_diri, f'*{datei}*'))
    # Filter files containing grouped name in their names
    datei_Groupedfilelist = [filename for filename in datei_filelist if GroupedName in filename]
    # Sort the filelist using the custom sorting key
    sorted_datei_filelist = sorted(datei_Groupedfilelist, key=extract_numeric_part)
    
    if len(sorted_datei_filelist)==0:
        print('No data for',datei,GroupedName)
    else:
        # proceed
        alldates_files = alldates_files+sorted_datei_filelist

#=====1. For each date, merge across all vertical layers======
#--------------------------------------------------------------------------

# Read in the grouped names for all times given
for datei in GivenDays_fmt2:
    # Get a list of file paths in the specified directories with the given date
    datei_filelist = glob.glob(os.path.join(h2files_diri, f'*{datei}*'))
    # Filter files containing grouped name in their names
    datei_Groupedfilelist = [filename for filename in datei_filelist if GroupedName in filename]
    # Sort the filelist using the custom sorting key
    sorted_datei_filelist = sorted(datei_Groupedfilelist, key=extract_numeric_part)
    if len(sorted_datei_filelist)==0:
        print('No data for',datei,GroupedName)
    else:
        # proceed to merge for each date, all vertical layers
        # use the first variable to test if datei has been processed already
        Testoutput_fileName = casename+'.cam.h2.MassConserve_latlon015.'+var_dic[GroupedName][0]+'.'+datei+'.nc'
        
        if os.path.exists(singleDay_diri+Testoutput_fileName):
            print(f"{datei}; has been processed.")
        else:
            print(f"Start to merge all vertical layers for: {datei}.")
            # get correct time array
            from datetime import datetime, timedelta
            datetime_list = []
            current_date = datetime(int(datei[:4]),int(datei[5:7]),int(datei[8:10]), 0, 0, 0)

            # Create a list of datetime values with a one-hour interval
            for hour in range(16, 23):
                dt = current_date.replace(hour=hour, minute=0, second=0, microsecond=0)
                datetime_list.append(dt)
            # Convert the datetime objects to their string representations
            formatted_datetime_list = [dt.strftime('%Y-%m-%dT%H:%M:%S.%f') for dt in datetime_list]

            # Loop through the list of file paths and read each variable into a separate dataset
            for fileidx in range(len(sorted_datei_filelist)):
                file_path = sorted_datei_filelist[fileidx]
                # print(file_path) # if to check which files went into the loop
                
                dataset = xr.open_dataset(file_path)
                levidx = file_path.split('.')[-2]

                # Add 'lev' as a new coordinate
                dataset = dataset.assign_coords(lev=levidx)
                # Expand all variables to include the 'lev' dimension
                for var_name in var_dic[GroupedName]:
                    dataset[var_name] = dataset[var_name].expand_dims('lev', axis=-1)

                # Merge the datasets along the 'lev' dimension
                if fileidx==0:
                    combined_dataset = dataset
                else:
                    # # Add each variable from the dataset to the combined dataset
                    combined_dataset = xr.concat([combined_dataset, dataset], dim='lev')

            # after merging
            # Correct for the time array - Select times around TROPOMI overpass T16:00:00 to T22:00:00
            combined_dataset = combined_dataset.assign_coords(time=formatted_datetime_list)
            # Make sure 'time' is in the correct format
            combined_dataset['time'] = combined_dataset['time'].astype('datetime64[ns]')

            # # save the combined 
            # combined_dataset.to_netcdf(singleDay_diri+output_fileName)
            # print('Saved to:',singleDay_diri+output_fileName)
            
            # each variable write into a .nc file
            for varname in var_dic[GroupedName]:
                var_ds = combined_dataset[[varname]]
                output_fileName = casename+'.cam.h2.MassConserve_latlon015.'+varname+'.'+datei+'.nc'
                var_ds.to_netcdf(singleDay_diri+output_fileName)
                print('Saved to:',singleDay_diri+output_fileName)
            
#=====2. If then merge across all times for each variable======
#--------------------------------------------------------------------------
# Merge for Time for each vertical layered variables
for varname in var_dic[GroupedName]:
    # Check and get file names
    output_fileName = casename+'.cam.h2.MassConserve_latlon015.'+varname+'.'+GivenDays[0]+'T'+GivenDays[-1]+'.nc'
    
    if os.path.exists(FileOutByVar_dir+output_fileName):
        print(f"{varname}; has been saved.")
    else:
        alldates_files = []
        # Read in the grouped names for all times given
        for datei in GivenDays_fmt2:
            # Get a list of file paths in the specified directories with the given date
            datei_filelist = glob.glob(os.path.join(MergedVerticalh2files_diri, f'*{datei}*'))
            # Filter files containing grouped name in their names
            varstr =f'.{varname}.' # add the dot so the selection is effective on filename
            datei_Groupedfilelist = [filename for filename in datei_filelist if varstr in filename]
            # Sort the filelist using the custom sorting key
            sorted_datei_filelist = sorted(datei_Groupedfilelist, key=extract_numeric_part)

            if len(sorted_datei_filelist)==0:
                print('No data for',datei,varname)
            else:
                # proceed
                alldates_files = alldates_files+sorted_datei_filelist

        print(len(alldates_files))

        ### Have to merge across times and vertical layers
        # simply read in all together
        combined_dataset = xr.open_mfdataset(alldates_files)
        # Write for all given dates into a single .nc file
        combined_dataset.to_netcdf(FileOutByVar_dir+output_fileName)
        print('Saved to:',FileOutByVar_dir+output_fileName)
