'''
This code is designed to merge all one-vertical-layer variables ('NoVertvars') for each day for re-gridded SE outputs of MUSICA-V0 in regular 0.15x0.15 lat and lon grids using mass conservative method

The MUSICAv0 output has already been mass-conservatively regridded into a resolution of 0.15x0.15 lat-lon but in separate files for each date

MODIFICATION HISTORY:
    Madankui Tao, 17, December, 2023: VERSION 1.0
    - Initial version
'''

GroupedName = 'NoVertvars'
# single species files use .XX. to make sure the file selection was effective

startDate = "2018-07-01"
endDate = "2018-08-01" 
endDateExclude = "2018-08-02" # the created datearray exclude this date

var_dic = {
           'NoVertvars':['TROP_P','PS'],
          }

# casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1base01'
# casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1basehourlyNEI2017'
# casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1basehourlyNOotherJulyMeanNEI2017'
# casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1baseJulyMeanNEI2017'
### with 6-hr nudging
casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.6HrNudgeTS1hourlyNEI2017'

# directory specification
h2files_diri = '/net/fs09/d0/taoma528/CESM22/Regridded_MUSICA_Output/2018_1330LT_TROPOMIcomp/MassConserve_latlon015_MUSICAoutput/'+casename+'/h2/'
singleDay_diri = '/net/fs09/d0/taoma528/CESM22/Regridded_MUSICA_Output/2018_1330LT_TROPOMIcomp/MassConserve_latlon015_MUSICAoutput/'+casename+'/h2_VertiVarsByVarDate/'
FileOutByVar_dir = '/net/fs09/d0/taoma528/CESM22/Regridded_MUSICA_Output/2018_1330LT_TROPOMIcomp/MassConserve_latlon015_MUSICAoutput/'+casename+'/h2_ByVar/'
MergedVerticalh2files_diri = '/net/fs09/d0/taoma528/CESM22/Regridded_MUSICA_Output/2018_1330LT_TROPOMIcomp/MassConserve_latlon015_MUSICAoutput/'+casename+'/h2_VertiVarsByVarDate/'

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
# print(alldates_files)


#--------------------------------------------------------------------------
# Fix time coordinates
from datetime import datetime, timedelta
startDate = datetime(int(startDate[:4]),int(startDate[5:7]),int(startDate[8:10]))
endDate = datetime(int(endDate[:4]),int(endDate[5:7]),int(endDate[8:10]))

# Initialize an empty list to store all datetime values
datetime_list = []

# Iterate through the date range
current_date = startDate
while current_date <= endDate:
    # For each date, generate times from T16 to T22
    for hour in range(16, 23):
        dt = current_date.replace(hour=hour, minute=0, second=0, microsecond=0)
        datetime_list.append(dt)
    current_date += timedelta(days=1)

# Convert the datetime objects to their string representations
formatted_datetime_list = [dt.strftime('%Y-%m-%dT%H:%M:%S.%f') for dt in datetime_list]
# print(formatted_datetime_list[0],formatted_datetime_list[-1])

#=====Merge across all times for each variable======
#--------------------------------------------------------------------------
### Only need to merge across time
if GroupedName=='NoVertvars':
    # simply read in all together
    combined_dataset = xr.open_mfdataset(alldates_files)
    # Re-write the time dimension to the correct value
    # Select times around TROPOMI overpass T16:00:00 to T22:00:00
    combined_dataset = combined_dataset.assign_coords(time=formatted_datetime_list)
    # Make sure 'time' is in the correct format
    combined_dataset['time'] = combined_dataset['time'].astype('datetime64[ns]')
    
    # each variable write into a .nc file
    for varname in var_dic[GroupedName]:
        var_ds = combined_dataset[[varname]]
        output_fileName = casename+'.cam.h2.MassConserve_latlon015.'+varname+'.'+GivenDays[0]+'T'+GivenDays[-1]+'.nc'
        var_ds.to_netcdf(FileOutByVar_dir+output_fileName)
        print('Saved to:',FileOutByVar_dir+output_fileName)