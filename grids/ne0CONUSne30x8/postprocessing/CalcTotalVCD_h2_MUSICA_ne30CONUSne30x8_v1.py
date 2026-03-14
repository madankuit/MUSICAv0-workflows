# This script is used for CESM2.2 CAM-chem with SE dynamical core simulations 
'''
this code is designed to calculate total vertical column densities of HCHO or NO2 or CO of h2 MUSICA-V0 outputs in model output ne30CONUSne30x8 grids over the CONUS

process separately for each variable, for the duration of the given period

MODIFICATION HISTORY:
    Madankui Tao, May, 29, 2024: VERSION 1.0
    - Based on CalcTropVCD_h2_MUSICA_ne30CONUSne30x8_v2.py 
'''
#================================================================================================
# inputs
varname = 'CO' # 'HCHO' or 'NO2'
# Variables to modify
startdate = "2018-07-01" #"2018-06-30"
enddate = "2018-07-05" # inclusive of last day

# varname = 'HCHO' # 'HCHO' or 'NO2'
# # Variables to modify
# startdate = "2018-07-01" #"2018-06-30"
# enddate = "2018-08-01" # inclusive of last day

FM1_startdate = startdate.replace("-", "")
FM1_enddate = enddate.replace("-", "")

# ### One-month NEI simulations
# casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1base01'
# casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1baseJulyMeanNEI2017'
# casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1basehourlyNEI2017'
# casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1basehourlyNOotherJulyMeanNEI2017'
# casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1basedailyNOotherJulyMeanNEI2017'

### 5-day sensitivity simulations
# casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.July5DaySensiTest.TS1base01'
# casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.July5DaySensiTest.TS1Anthro70Perct01'
casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.July5DaySensiTest.TS1BioEmis70Perct01'

# where to put the output files
MUSICAVCDSave_diri = f'/net/fs09/d0/taoma528/CESM22/Calculated_MUSICA_VCD/TotalVCD/{casename}/'
SavePath = f'{MUSICAVCDSave_diri}{casename}.cam.h2.TotalVCD.{varname}.{FM1_startdate}T{FM1_enddate}.nc'

MetPath = '/net/fs09/d0/taoma528/CESM22/Calculated_MUSICA_VCD/TroposphericVCD/met_f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1base01.nc'

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
import sys
sys.path.append('/home/taoma528/Scripts/CESM_analysis/functions')
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
svante_archive = '/net/fs09/d0/taoma528/CESM22/archive/'
RunPath = svante_archive+casename+'/atm/hist/'

# for hourly average file
histfreq = 'h2'

# need to save CO, HCHO, NO2, and also the meteorological variables and tropopause height (to calculate VCD)
variables_to_select = ['CO','CH2O', 'NO2','lat','lon','date','area','hyam','hybm',]#'T','TROP_P']#,'PS','PMID']

# Get a list of file paths in the specified directories with the given keyword; all dates
Run_filelist = glob.glob(os.path.join(RunPath, f'*{histfreq}*'))
Run_filelist = sorted(Run_filelist)

### Read one file in as an example
testMUSICAv0_ds = xr.open_dataset( Run_filelist[0] )[variables_to_select]
numlevs = len(testMUSICAv0_ds.lev)
ncols = len(testMUSICAv0_ds.ncol.values)
    
lat_da = testMUSICAv0_ds.lat
lon_da = testMUSICAv0_ds.lon

met_ds = xr.open_dataset(MetPath) # as 'PS','PMID' are not saved

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
        elif varname in ['NO2','CO']:
            houri_CONUS_var_da = houri_MUSICAv0_ds[varname].where(mask)
        
        # try to get the index of CONUS pixels | horizontal
        example_ar = houri_CONUS_var_da.isel(lev=0).values
        # Find indices of non-NaN values
        CONUS_indices = np.where(~np.isnan(example_ar))[0]
        
        ### Loop through pixels
        # Loop through each ncol (pixel) | Only calculate for the CONUS pixels
        for colidx in CONUS_indices:
            # -----------------------------
            # TOTAL column (all levels; PMID midpoints; TOA->surface order)
            # -----------------------------
            SurfP_Pa = float(houri_met_ds["PS"].sel(ncol=colidx).values)          # Pa
            PMID_all_Pa = houri_met_ds["PMID"].sel(ncol=colidx).values            # (lev,) Pa midpoints, TOA->sfc
            T_all_K     = houri_met_ds["T"].sel(ncol=colidx).values               # (lev,) K, TOA->sfc
            var_all_molmol = houri_CONUS_var_da.sel(ncol=colidx).values           # (lev,) mol/mol, TOA->sfc

            # basic finite filter only
            valid = np.isfinite(PMID_all_Pa) & np.isfinite(T_all_K) & np.isfinite(var_all_molmol)
            PMID_ncoli_Pa = PMID_all_Pa[valid]
            T_ncoli_K     = T_all_K[valid]
            var_ncoli_molmol = var_all_molmol[valid]

            nlayer = len(PMID_ncoli_Pa)

            # if something went wrong (shouldn't), skip
            if nlayer < 2:
                Datei_hourly_VCD_ar[houridx, colidx] = np.nan
            else:
                Avog = 6.0221409e23   # molecules / mol
                R = 8.314462          # Pa m^3 mol^-1 K^-1
                DryAirMass = 0.029    # kg/mol
                m3Tocm3 = 1e6

                ncoli_var_molec_cm2 = []

                # Loop from surface -> TOA (reverse the TOA->surface arrays)
                for ridx in range(nlayer - 1, -1, -1):

                    p_i = PMID_ncoli_Pa[ridx]          # Pa
                    T_i = T_ncoli_K[ridx]              # K
                    x_i = var_ncoli_molmol[ridx]       # mol/mol

                    # Thickness (cm) using your Pres_TO_Height(P_mid, P_sfc)
                    # - bottom layer (closest to surface): from midpoint to surface
                    # - other layers: difference between adjacent midpoints
                    # NOTE: top-of-atmos above top midpoint is ignored (as you intended).
                    if ridx == (nlayer - 1):
                        dz_cm = Pres_TO_Height(p_i, SurfP_Pa) * 1e2
                    else:
                        p_below = PMID_ncoli_Pa[ridx + 1]
                        dz_cm = (Pres_TO_Height(p_i, SurfP_Pa) - Pres_TO_Height(p_below, SurfP_Pa)) * 1e2

                    # Air density (kg/m^3)
                    rho_air = (p_i * DryAirMass) / (R * T_i)

                    # molecules/cm^2 in layer
                    vcd_layer = x_i * Avog * (rho_air / DryAirMass) * (dz_cm / m3Tocm3)

                    # add to the array to calculate the sum of all levels
                    ncoli_var_molec_cm2.append(vcd_layer)

                Datei_hourly_VCD_ar[houridx, colidx] = np.nansum(np.asarray(ncoli_var_molec_cm2))

    ### Add to AllDates array
    AllDates_hourly_VCD_ar[Dateidx*24:(Dateidx+1)*24,:] = Datei_hourly_VCD_ar
    print(f'Finished calculation for {Datei}')
    
#================================================================================================
## Write to Xr Dataset and save
import datetime
today_date = datetime.date.today()

# write this variable to a dataarray
MUSICA_TotalVCD_ds = xr.Dataset(
                    {
                        varname+'_TotalVCD': (['time','ncol'], AllDates_hourly_VCD_ar),
                    },
                    coords={'time': ('time', AllDates_UTC_hours_ar), 
                            'ncol': ('ncol', MUSICAv0_ds.ncol.values), 
                            },
                    attrs={'process date': str(today_date.strftime("%Y-%m-%d")),
                           'description': 'MUSICA total VCD over the CONUS',
                          },
                    )

# Add the attributes of coordinates
for var_name in ['time','ncol']:    
    # Copy attributes from source to target variable
    for attr_name, attr_value in MUSICAv0_ds[var_name].attrs.items():
        MUSICA_TotalVCD_ds[var_name].attrs[attr_name] = attr_value

MUSICA_TotalVCD_ds[varname+'_TotalVCD'].attrs['longname'] = 'total vertical column density of '+varname
MUSICA_TotalVCD_ds[varname+'_TotalVCD'].attrs['unit'] = 'molec/cm2'

### Add lat and lon DataArray
MUSICA_TotalVCD_ds = MUSICA_TotalVCD_ds.assign(lat=lat_da)
MUSICA_TotalVCD_ds = MUSICA_TotalVCD_ds.assign(lon=lon_da)

#================================================================================================
# save to .nc
MUSICA_TotalVCD_ds.to_netcdf(SavePath)
print('Save to:',SavePath)

        
        

