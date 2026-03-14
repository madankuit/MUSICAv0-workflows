# This script is used for CESM2.2 CAM-chem with SE dynamical core simulations
'''
this code is designed to regrid h2 (hourly) SE outputs of MUSICA-V0 to regular 0.15x0.15 lat and lon grids using mass conservative method

Based on 'ConservativeRegrid_MUSICAOutputs.ipynb' from MUSICA Tutorial Nov. 2021

Regrid only the horizontal resolution to regular lat,lon 

MODIFICATION HISTORY:
    Madankui Tao, 4, December, 2023: VERSION 1.0
    - Initial version
    Madankui Tao, 16, December, 2023: VERSION 1.1
    - Add CO to verticalvars, adjust for different casename
'''
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

# Functions copied from Cheyenne
import sys
sys.path.insert(0,'/home/taoma528/Scripts/CESM_analysis/MergeNEI2017_to_CAMS/NCAR_packages/')
from dsj.array.chk import chk
from dsj.plot.Plot_2D import Plot_2D
from dsj.time.expand_date import Expand_date
from dsj.analysis.Regridding_ESMF import Add_bounds, Regridding

#================================================================================================
### Providing all files needed to regrid
# ne0CONUSne30x8 grid file as the source
SERR_scrip_file = '/net/fs09/d0/taoma528/CESM22/grids/ne0CONUS_ne30x8_np4_SCRIP.nc'

# Target Grid information file for the FV grid
ESMFmap_015grid_file = '/net/fs09/d0/taoma528/CESM22/grids/FV_gridinfo_0.15_c20231204.nc'
# ESMFmap_grid_file = '/net/fs09/d0/taoma528/CESM22/grids/CheyCopy_FV_gridinfo_015_c20231202.nc'

# Weights file
Regridding_015weights_file = '/net/fs09/d0/taoma528/CESM22/grids/ne0CONUSne30x8_ESMFmap_0.15x0.15_cubit_conserve_cams_c20231204.nc'

#================================================================================================
# Read in MUSICA outputs
# specify the MUSICA directory
svante_archive = '/net/fs09/d0/taoma528/CESM22/archive/'

# casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1base01'
# casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1basehourlyNEI2017'
# casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1basehourlyNOotherJulyMeanNEI2017'
# casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1baseJulyMeanNEI2017'
### with 6-hr nudging
casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.6HrNudgeTS1hourlyNEI2017'

RunPath = svante_archive+casename+'/atm/hist/'

# for hourly average file
histfreq = 'h2'
# where to put the regridded files
histfreqfiles_diri = '/net/fs09/d0/taoma528/CESM22/Regridded_MUSICA_Output/2018_1330LT_TROPOMIcomp/MassConserve_latlon015_MUSICAoutput/'+casename+'/'+histfreq+'/'

# need to save HCHO, NO2, and also the meteorological variables and tropopause height (to calculate VCD)
# variables_to_select = ['CO','CH2O', 'NO2','lat','lon','date','area','hyam','hybm','T','TROP_P','PS','PMID']

NoVertvars = ['TROP_P','PS']
vertivars = ['CO','CH2O', 'NO2','T','PMID']

# Get a list of file paths in the specified directories with the given keyword
Run_filelist = glob.glob(os.path.join(RunPath, f'*{histfreq}*'))
Run_filelist = sorted(Run_filelist)

#================================================================================================
# Regrid for variables that do NOT have vertical information

# Uncomment if test with one file
# fileidxs = [0]
# for fileidx in fileidxs:

# otherwise use range(len(Run_filelist))
for fileidx in range(len(Run_filelist)):
    Run_filepathi = Run_filelist[fileidx]
    Run_filenamei = Run_filepathi.split('/')[-1]
    
    # read in hourly data
    MUSICAv0_ds = xr.open_dataset( Run_filepathi )
    # Select times around TROPOMI overpass T16:00:00 to T22:00:00
    selTimeMUSICAv0_ds = MUSICAv0_ds.isel(time=slice(15,22))
    
    numlevs = len(selTimeMUSICAv0_ds.lev)

    ### Write Variables with no Vertical levels together
    print( 'Regridding: ', Run_filenamei, "['TROP_P','PS']" ) #, datetime.now()
    print( '************************************************************************' )
    dst_file = histfreqfiles_diri+Run_filenamei[:-3]+'.MassConserve_latlon015.NoVertvars.nc'
    if os.path.exists(dst_file):
        print(f"The file {dst_file} exists.")
    else:
        rr = Regridding( selTimeMUSICAv0_ds[NoVertvars], src_grid_file=SERR_scrip_file, dst_grid_file=ESMFmap_015grid_file, 
                         wgt_file=Regridding_015weights_file, method='Conserve', fields=NoVertvars,
                    dst_file=dst_file, save_wgt_file=False, save_results=True, check_results=False, 
                         check_timings=True, creation_date=False, nc_file_format='NETCDF4' )
    print( 'Write to: ', dst_file )
    
    ### Regrid for variables that HAVE vertical information
    print( 'Regridding: ', Run_filenamei, "['CO','CH2O', 'NO2','T','PMID']" ) #, datetime.now()
    print( '************************************************************************' )
    # Each vertical lev into one file to combine and re-write later
    for levidx in range(numlevs):
        dst_file = histfreqfiles_diri+Run_filenamei[:-3]+'.MassConserve_latlon015.vertivars.'+str(levidx)+'.nc'
        if os.path.exists(dst_file):
            print(f"The file {dst_file} exists.")
        else:
            rr = Regridding( selTimeMUSICAv0_ds[vertivars].isel(lev=levidx), src_grid_file=SERR_scrip_file, dst_grid_file=ESMFmap_015grid_file, 
                         wgt_file=Regridding_015weights_file, method='Conserve', fields=vertivars,
                    dst_file=dst_file, save_wgt_file=False, save_results=True, check_results=False, 
                         check_timings=True, creation_date=False, nc_file_format='NETCDF4' )
    print( 'Saved all layers' )

