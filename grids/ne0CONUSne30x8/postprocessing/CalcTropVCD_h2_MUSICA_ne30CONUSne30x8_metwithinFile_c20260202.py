# This script is used for CESM2.2 CAM-chem with SE dynamical core simulations 
'''
this code is designed to calculate tropospheric vertical column densities of HCHO or NO2 of h2 MUSICA-V0 outputs in model output ne30CONUSne30x8 grids over the CONUS

process separately for HCHO and NO2, for the duration of the given period

Can run with sbatch run_CalcTropVCD.slurm to have it running in fdr in the background

MODIFICATION HISTORY:
    M. Tao, Feb, 2, 2026: VERSION 1.1
    - Copied from ACP_MUSICANEI_scripts/CalcTropVCD_h2_MUSICA_ne30CONUSne30x8_v2.py
    - Modified to use meteorology variables saved in the h2 file, instead of read externally
'''
# ============================================================
# USER CONFIGURATION — set these paths before running
# ============================================================
# Path to your CESM archive directory, e.g. '/path/to/CESM/archive/'
svante_archive = ''

# Path to directory for saving tropospheric VCD output files
# e.g. '/path/to/Calculated_MUSICA_VCD/TroposphericVCD/<casename>/'
MUSICATropVCDSave_diri = ''

# casename options
casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.BASE.From20230801.01'
# casename = ''
# ============================================================

#================================================================================================
# inputs
varname = 'HCHO' # 'HCHO' or 'NO2'
# Variables to modify | calculation is so slow, have to break into ~ 1month/calculation
startdate = "2024-06-01" #
enddate = "2024-06-30" # inclusive of last day
# enddate = "2024-09-01" # inclusive of last day

FM1_startdate = startdate.replace("-", "")
FM1_enddate = enddate.replace("-", "")

SavePath = f'{MUSICATropVCDSave_diri}{casename}.cam.h2.TroposphericVCD.{varname}.{FM1_startdate}T{FM1_enddate}.nc'

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

FMT1_days = []
MUSICA_days = []
import datetime
Read_start = datetime.datetime.strptime(startdate, "%Y-%m-%d")
Read_end_plus1 = datetime.datetime.strptime(enddate, "%Y-%m-%d") + timedelta(days=1)
date_generated = [Read_start + timedelta(days=x) for x in range(0, (Read_end_plus1-Read_start).days)]

# In the format of filename date
for date in date_generated:
    FMT1_days.append(date.strftime("%Y%m%d"))
    MUSICA_days.append(date.strftime("%Y-%m-%d"))
    
# in a dictionary to adjust for different date format
FMT1date_onMUSICAdate = dict(zip(FMT1_days,MUSICA_days))

#================================================================================================
### Check if the file already exists
if os.path.exists(SavePath):
    print(f"The file {SavePath} exists.")
    # exit the script
    sys.exit()
else:
    print(f'Start to process for {varname} to store to {SavePath}')
    
#================================================================================================
# For a given case: Read hourly model outputs h2 files
# specify the MUSICA directory
RunPath = svante_archive+casename+'/atm/hist/'

# for hourly average file
histfreq = 'h2'

# need to save HCHO, NO2, and also the meteorological variables and tropopause height (to calculate VCD)
variables_to_select = ['CH2O', 'NO2','lat','lon','date','area','hyam','hybm',]#'T','TROP_P']#,'PS','PMID']

# Get a list of file paths in the specified directories with the given keyword; all dates
Run_filelist = glob.glob(os.path.join(RunPath, f'*{histfreq}*'))
Run_filelist = sorted(Run_filelist)

### Read one file in as an example
testMUSICAv0_ds = xr.open_dataset( Run_filelist[0] )[variables_to_select]
numlevs = len(testMUSICAv0_ds.lev)
ncols = len(testMUSICAv0_ds.ncol.values)
    
lat_da = testMUSICAv0_ds.lat
lon_da = testMUSICAv0_ds.lon

met_vars = ['TROP_P','PS','PMID','T']

#================================================================================================
# Size of final dataset to store VCD (end product)
AllDates_hourly_VCD_ar = np.zeros([len(MUSICA_days)*24,ncols])
AllDates_hourly_VCD_ar[:,:] = np.nan

# store the hours
AllDates_UTC_hours_ar = np.empty((0,), dtype='datetime64[ns]')

#================================================================================================
#================================================================================================
### Loop through days
for Dateidx,Datei in enumerate(MUSICA_days):
    Run_filepathi = next((file for file in Run_filelist if Datei in file), None)
    Run_filenamei = Run_filepathi.split('/')[-1]
    
    # read in hourly data for the selected variables
    MUSICAv0_ds = xr.open_dataset( Run_filepathi )[variables_to_select]
    met_ds = xr.open_dataset( Run_filepathi )[met_vars]
    
    # numlevs = len(MUSICAv0_ds.lev)
    # ncols = len(MUSICAv0_ds.ncol.values)

    # Use the base case to get the adjusted lons and add to each case
    #-----------------------------------------------------------------------
    # MUSICA lons goes from 0-360, convert it to +-180 | Need to adjust lon_right and lon_left for MUSICA!!
    mdllon = MUSICAv0_ds.isel(time=1)['lon'].values
    mdlats = MUSICAv0_ds.isel(time=1)['lat'].values
    Adjustedlons = np.where(mdllon >= 180, mdllon - 360, mdllon)
    # Create a DataArray for the adjusted lon
    Adjustedlons_da = xr.DataArray(
                                    Adjustedlons,  
                                    dims=('ncol'),
                                    coords={'ncol': MUSICAv0_ds.ncol.values} 
                                    )

    # Add the DataArray back to the dataset
    MUSICAv0_ds['Adjustedlons'] = Adjustedlons_da

    ### Only do calculations over the CONUS domain
    # extract for a given region
    fileregion = 'CONUS' 
    lon_right,lat_bot,lon_left,lat_up = latlonbound(fileregion)

    # Use the base case to create a mask based on latitude and longitude
    mask = (MUSICAv0_ds['Adjustedlons'] >= lon_right) & \
           (MUSICAv0_ds['Adjustedlons'] <= lon_left) & \
            (MUSICAv0_ds['lat'] >= lat_bot) & \
           (MUSICAv0_ds['lat'] <= lat_up) 
    
    # add the times
    Datei_hours = MUSICAv0_ds.time.values
    AllDates_UTC_hours_ar = np.append(AllDates_UTC_hours_ar, Datei_hours)
    
    #-----------------------------------------------------------------------
    # loop through hours and calculate column densities, store to AllDates_hourly_VCD_ar[Dateidx*24:(Dateidx+1)*24,:]
    Datei_hourly_VCD_ar = np.zeros([len(Datei_hours),ncols])
    Datei_hourly_VCD_ar[:,:] = np.nan
    
    ### Loop through each hour
    for houridx, houri in enumerate(Datei_hours):
        houri_MUSICAv0_ds = MUSICAv0_ds.sel(time = houri)
        houri_met_ds = met_ds.sel(time = houri)
        # mask for only CONUS
        if varname=='HCHO':
            houri_CONUS_var_da = houri_MUSICAv0_ds['CH2O'].where(mask)
        elif varname=='NO2':
            houri_CONUS_var_da = houri_MUSICAv0_ds[varname].where(mask)
        
        # try to get the index of CONUS pixels | horizontal
        example_ar = houri_CONUS_var_da.isel(lev=0).values
        # Find indices of non-NaN values
        CONUS_indices = np.where(~np.isnan(example_ar))[0]
        
        ### Loop through pixels
        # Loop through each ncol (pixel) | Only calculate for the CONUS pixels
        for colidx in CONUS_indices:
            # select conditions at this pixel
            tropopauseP_Pa = float(houri_met_ds['TROP_P'].sel(ncol=colidx).values) # tropopause pressure in Pa
            SurfP_Pa = float(houri_met_ds['PS'].sel(ncol=colidx).values) # surface pressure in Pa
            allPMID_Pa = houri_met_ds['PMID'].sel(ncol=colidx) # vertical pressure in Pa
            temp_K = houri_met_ds['T'].sel(ncol=colidx) # temperature in K

            # get the index of lev of PMID that is lower than tropopause (Pressure > Tropopause pressure)
            # Find the indices where allPMID_Pa is greater than tropopauseP_Pa
            verticalTropindices = np.where(allPMID_Pa.values >= tropopauseP_Pa)[0] # to an array

            ### Get everything below the tropopause (ordered from tropopause to the surface)
            PMID_ncoli_Pa = allPMID_Pa.isel(lev=slice(verticalTropindices[0],verticalTropindices[-1]+1)).values
            T_ncoli_K = temp_K.isel(lev=slice(verticalTropindices[0],verticalTropindices[-1]+1)).values

            # variable value
            varval_ncoli_molPmol = houri_CONUS_var_da.sel(ncol=colidx).isel(lev=slice(verticalTropindices[0],verticalTropindices[-1]+1)).values

            # number of layers
            troplayernum = len(verticalTropindices) 

            #####
            # # calculate for the given variable, output sum up
            ncoli_var_molec_cm2 = []

            ## Loop through the vertical layers within the troposphere, 
            for LAYidx, reversedLAYidx in enumerate(range(troplayernum - 1, -1, -1)):
                # backwards from surface to Tropopause (currently everything is order from TOA to surface)
                MUSICA_levi_PMID = PMID_ncoli_Pa[reversedLAYidx] # Pa
                MUSICA_levi_TEMP = T_ncoli_K[reversedLAYidx] # K
                var_moleFraci = varval_ncoli_molPmol[reversedLAYidx] # mole var/mole air
                # print(MUSICA_levi_PMID,MUSICA_levi_TEMP,var_moleFraci)

                if reversedLAYidx==(troplayernum-1):
                    # If for the near-surface level, use the height to the surface
                    MUSICA_levi_H_cm = Pres_TO_Height(MUSICA_levi_PMID,SurfP_Pa)*(1e2)
                else:
                    # if not at the near-surface level, get the pressure of the layer below as well and calculate the diff
                    MUSICA_levibelow_PMID = PMID_ncoli_Pa[reversedLAYidx+1] # Pa float(MUSICA_levibelow_ds.PMID.values)
                    MUSICA_levi_H_cm = Pres_TO_Height(MUSICA_levi_PMID,SurfP_Pa)*(1e2)-Pres_TO_Height(MUSICA_levibelow_PMID,SurfP_Pa)*(1e2)

                # Calculate the sum of all cells across a layer (unit originally in mol/mol)
                Avog = 6.0221409e23 # molec/mole
                m3Tocm3 = 1e6
                R = 8.314462 # kg⋅m2⋅s−2⋅K−1⋅mol−1 / m3⋅Pa⋅K−1⋅mol−1
                DryAirMass = 0.029 # kg/mol
                rho_airi = (MUSICA_levi_PMID*DryAirMass)/(R*MUSICA_levi_TEMP)

                levi_var_molec_cm2 = var_moleFraci*Avog*(1/DryAirMass)*rho_airi*(1/m3Tocm3)*MUSICA_levi_H_cm

                # add to the array to calculate the sum of all levels
                ncoli_var_molec_cm2.append(levi_var_molec_cm2)

            # add tropVCD summed up for all layers, put to the corresponding ncol grid
            Datei_hourly_VCD_ar[houridx,colidx] = np.nansum(np.array(ncoli_var_molec_cm2))
        
    ### Add to AllDates array
    AllDates_hourly_VCD_ar[Dateidx*24:(Dateidx+1)*24,:] = Datei_hourly_VCD_ar
    print(f'Finished calculation for {Datei}')
    
#================================================================================================
## Write to Xr Dataset and save
import datetime
today_date = datetime.date.today()

# write this variable to a dataarray
MUSICA_TropVCD_ds = xr.Dataset(
                    {
                        varname+'_TropVCD': (['time','ncol'], AllDates_hourly_VCD_ar),
                    },
                    coords={'time': ('time', AllDates_UTC_hours_ar), 
                            'ncol': ('ncol', MUSICAv0_ds.ncol.values), 
                            },
                    attrs={'process date': str(today_date.strftime("%Y-%m-%d")),
                           'description': 'MUSICA tropospheric VCD over the CONUS',
                          },
                    )

# Add the attributes of coordinates
for var_name in ['time','ncol']:    
    # Copy attributes from source to target variable
    for attr_name, attr_value in MUSICAv0_ds[var_name].attrs.items():
        MUSICA_TropVCD_ds[var_name].attrs[attr_name] = attr_value

MUSICA_TropVCD_ds[varname+'_TropVCD'].attrs['longname'] = 'tropospheric vertical column density of '+varname
MUSICA_TropVCD_ds[varname+'_TropVCD'].attrs['unit'] = 'molec/cm2'

### Add lat and lon DataArray
MUSICA_TropVCD_ds = MUSICA_TropVCD_ds.assign(lat=lat_da)
MUSICA_TropVCD_ds = MUSICA_TropVCD_ds.assign(lon=lon_da)

#================================================================================================
# save to .nc
MUSICA_TropVCD_ds.to_netcdf(SavePath)
print('Save to:',SavePath)

        
        

