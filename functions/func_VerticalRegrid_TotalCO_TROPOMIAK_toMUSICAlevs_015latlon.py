'''
This function is used to re-grid TROPOMI AK already processed to 0.15x0.15 degree to the vertical resolution of MUSICA for a given region on datei (YYYYMMDD)

# Revised from func_VerticalRegrid_TROPOMIAK_toMUSICAlevs.py, modified given different ways of processing MUSICA and TROPOMI
# Note that TROPOMI CO is also recorded from TOA to surface based on layers, same as MUSICA

MODIFICATION HISTORY:
    M. Tao, 13, December, 2023: VERSION 1.0
    - Initial version
    M. Tao, 16, December, 2023: VERSION 1.1
    - Adjust for different casename
    M. Tao, 23, Feb, 2023: VERSION 2.0
    - Adjust for TotalCO (total VCD contains more layers than trop VCD)
'''

#=====Dependent library & functions======
#--------------------------------------------------------------------------
import numpy as np # for array manipulation and basic scientific calculation
import xarray as xr # To read NetCDF files
import pandas as pd

import warnings
# Ignore a specific warning
warnings.filterwarnings('ignore', message='Some warning message')
# Ignore the FutureWarning caused by iteritems()
warnings.filterwarnings("ignore", category=FutureWarning, message=".*iteritems.*")

import os
import glob

# interpolation
from scipy.interpolate import griddata # Simple regridding
from netCDF4 import Dataset # To write NetCDF files

### Calculate the mean concentrations of each species averaged at each region
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from func_MUSICA_DefineRegion import *

#=====Define the locations of TRROPOMI data======
#--------------------------------------------------------------------------
# Dictionaries for each species ['HCHO','NO2','CO']
# file header 
fileheader_dic = {'HCHO':'S5P_RPRO_L2__HCHO___',
                   'NO2':'S5P_RPRO_L2__NO2____',
                   'CO':'S5P_RPRO_L2__CO_____',
                  }

# Where to store the processed files
# Set your paths to the regridded TROPOMI data directories, e.g.:
L3_015TROPOMI_diri_dic = {'HCHO':'/path/to/TROPOMI/data/S5P_L2__HCHO___HiR_2/regrid_2018_MUSICA015/',
                           'NO2':'/path/to/TROPOMI/data/S5P_L2__NO2____HiR_2/regrid_2018_MUSICA015/',
                           'CO':'/path/to/TROPOMI/data/S5P_L2__CO_____HiR_2/regrid_2018_MUSICA015/',
                          }

#=====Example use======
#--------------------------------------------------------------------------
# # Example use
# casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1base01'
# datei = '20180701'
# fileregion = 'Northeast'

# vertically_gridded_AK_da = func_VerticalRegrid_TotalCO_TROPOMIAK_toMUSICAlevs_015latlon(casename, datei, fileregion)

#=====Dependent Function======
#--------------------------------------------------------------------------
def Height_TO_Pres(alt_m, P_surf):
    """This function calculates pressure (Pa) based on altitude (m) and surface pressure (P_surf)
        Scale height H defined as 8.5 km
    """
    H = 8.5  # scale height in km
    air_Pres = P_surf * np.exp(-alt_m / (H * 1e3))  # convert altitude to km
    return air_Pres

#=====Main Function======
#--------------------------------------------------------------------------

def func_VerticalRegrid_TotalCO_TROPOMIAK_toMUSICAlevs_015latlon(casename, datei, fileregion):
    """
    Calculate vertically gridded Average Kernel (AK) data to match MUSICA output.

    Parameters:
        casename (str): name of MUSICA case
        datei (str): The date to process the Average Kernel.
        fileregion (str): The region of interest.

    Returns:
        xarray.DataArray: The vertically gridded Average Kernel data for the given date and region.
    """
    
    varname = 'CO'
    #--------------------------------------------------------------------------
    # Get the target grid
    lon_right,lat_bot,lon_left,lat_up = latlonbound(fileregion)
    
    #--------------------------------------------------------------------------
    #--------------------------------------------------------------------------
    ### Use CO: also vertical layers from TOA (49500 m above) to surface (500 m above)
    # TROPOMI Test file for layers
    # TODO: set your path to the TROPOMI CO L2 directory, e.g. '/path/to/TROPOMI/S5P_L2__CO_____HiR_2/2018/'
    TROPOMI_L2COdiri = "/path/to/TROPOMI/S5P_L2__CO_____HiR_2/2018/"
    testCOfile = 'S5P_RPRO_L2__CO_____20180902T144647_20180902T162816_04600_03_020400_20220905T080034.nc'
    TROPOMI_CO_PRODUCT = xr.open_dataset(TROPOMI_L2COdiri+testCOfile,group="/PRODUCT/",engine="netcdf4")
    # TROP_CO_egda = xr.open_dataset(TROPOMI_L2COdiri+testCOfile,group="PRODUCT/SUPPORT_DATA/INPUT_DATA/").sel(time=0)
    layers_m = TROPOMI_CO_PRODUCT.layer.values
    layernum = len(layers_m)    
    # # get the array reversed
    # layers_m_reversed = np.flip(layers_m)
    #--------------------------------------------------------------------------
    # Calculated TM5 PMID for the given fileregion and datei, already processed to 0.15 degree lat-lon
    # TM5PMID_datei_layers have layers revserved from TOA to surface
    # TODO: set your path to the TM5 surface pressure data for CO, e.g.:
    # TM5SurfaceP_FileOut_diri = '/path/to/TM5MP_Model/'
    TM5SurfaceP_FileOut_diri = '/path/to/TM5MP_Model/'
    TM5SurfaceP_Path = TM5SurfaceP_FileOut_diri+'TM5MP_forCO_Model_SurfaceP_Global_015deg_2018-06-30T2018-08-02.nc'
    TM5SurfaceP_ds = xr.open_dataset(TM5SurfaceP_Path)
    regioniSurfaceP_ds = TM5SurfaceP_ds.sel(lat=slice(lat_bot,lat_up),lon=slice(lon_right,lon_left))
    dateiPs = regioniSurfaceP_ds.Ps.sel(time=datei).values
    # get the horizontal resolution grid
    regionilats = regioniSurfaceP_ds.lat
    regionilons = regioniSurfaceP_ds.lon 

    #--------------------------------------------------------------------------
    #--------------------------------------------------------------------------
    # Read in test MUSICA output file: get the vertical levels from MUSICA; does not matter which specific case
    # TODO: set your path to a representative MUSICA h2 output file, e.g.:
    # test_h2_MUSICAfile = '/path/to/MUSICA/archive/<casename>/atm/hist/<casename>.cam.h2.YYYY-MM-DD-SSSSS.nc'
    test_h2_MUSICAfile = '/path/to/MUSICA/archive/'+casename+'/atm/hist/'+casename+'.cam.h2.2018-07-01-03600.nc'
    testh2_ds = xr.open_dataset(test_h2_MUSICAfile)
    # hybrid level at midpoints: from TOA to Surface
    MUSICA_layers = testh2_ds.lev.values
    MUSICAlaynum = len(MUSICA_layers)
    del test_h2_MUSICAfile,testh2_ds
    #--------------------------------------------------------------------------
    # PMID processed to 0.15 degree already for the given case
    # TODO: set your path to regridded MUSICA output, e.g.:
    # regridByVar_diri = '/path/to/Regridded_MUSICA_Output/MassConserve_latlon015_MUSICAoutput/<casename>/h2_ByVar/'
    regridByVar_diri = '/path/to/Regridded_MUSICA_Output/MassConserve_latlon015_MUSICAoutput/'+casename+'/h2_ByVar/'
    MUSICA_PMID = xr.open_dataset( regridByVar_diri+casename+'.cam.h2.MassConserve_latlon015.PMID.20180701T20180801.nc' )
    regioni_MUSICA_PMID = MUSICA_PMID.PMID.sel(lat=slice(lat_bot, lat_up),lon=slice(lon_right,lon_left))
    # Approx at TROPOMI overpass time
    dateFMT2 = datei[:4]+'-'+datei[4:6]+'-'+datei[6:8]
    averageTimeStart,averageTimeEnd = timezone_approx1330LT(fileregion)
    sel_times = [np.datetime64(dateFMT2+averageTimeStart),np.datetime64(dateFMT2+averageTimeEnd)]
    # select for the approx two hours at 1 and 2 PM, then take the average
    regioni_MUSICA_PMID_1330LT = regioni_MUSICA_PMID.sel(time=sel_times).mean(dim='time') 

    #--------------------------------------------------------------------------
    #--------------------------------------------------------------------------
    # Calculate the absolute pressure in hPa | use the pressure at the middle of the interface level (match MUSICA) 
    TM5PMID_datei_layers = np.zeros([len(regionilats),len(regionilons),layernum])
    TM5PMID_datei_layers[:,:,:] = np.nan
    # calculate the pressure levels PMID using Ps_alldays_data from the last layer to use for the reversed 
    for LAYidx, reversedLAYidx in enumerate(range(layernum - 1, -1, -1)):
        # print(LAYidx,reversedLAYidx)    
        # Calculate the Pressure in Pa for each layer
        layi_m = layers_m[LAYidx]
        layi_Pa = Height_TO_Pres(layi_m, dateiPs)

        # put the revsered layer to the layer count from beginning
        TM5PMID_datei_layers[:,:,LAYidx] = layi_Pa

    # Save using reversed vertical layers, from TOA to Surface (Match MUSICA)
    TM5PMID_datei_ds = xr.DataArray(TM5PMID_datei_layers,
                                            dims=("lat", "lon","TOA_to_Surf_layer"),
                                            coords={
                                                "lat": regionilats,
                                                "lon": regionilons,
                                                "TOA_to_Surf_layer": layers_m,
                                            })

    #--------------------------------------------------------------------------
    # Read in TROPOMI AK for the given fileregion and datei; reverse the vertical layer order
    AK_FileOutPath = L3_015TROPOMI_diri_dic[varname]+fileheader_dic[varname]+'AK_Global_Regrid015deg_'+datei+'.nc'
    AK_ds = xr.open_dataset(AK_FileOutPath)
    # select for the date
    AK_xar = AK_ds['TROPOMI_'+varname+'_AK'].sel(time=datei)
    # CO AK should be the smallest near the surface    
    # select for the given region
    regioni_datei_AK_xar = AK_xar.sel(lat=slice(lat_bot, lat_up),lon=slice(lon_right,lon_left))    

    #--------------------------------------------------------------------------
    #--------------------------------------------------------------------------
    # Interpolate for each (lat,lon)
    vertically_gridded_AK = np.zeros([len(regionilats),len(regionilons),MUSICAlaynum])

    # reverse the order for the vertical layers
    for latidx in range(len(regionilats)):
        lati = regionilats[latidx]
        for lonidx in range(len(regionilons)):
            loni = regionilons[lonidx]
            # print(latidx,lonidx)
            # Values on the source grid
            values_source = regioni_datei_AK_xar.sel(lat=lati,lon=loni).values
            # Pressure levels on the source grid, used by TROPOMI in Pa
            pressure_source = TM5PMID_datei_ds.sel(lat=lati,lon=loni).values
            # Pressure levels on the target grid, used PMID by MUSICA in Pa
            pressure_target = regioni_MUSICA_PMID_1330LT.sel(lat=lati,lon=loni).values
            # Use numpy.interp to interpolate values from the source grid to the target grid
            values_target = np.interp(pressure_target, pressure_source, values_source)
            vertically_gridded_AK[latidx,lonidx,:] = values_target
            # print("Interpolated values on the target grid:", values_target)

    #--------------------------------------------------------------------------
    # Write to a xarray
    # Save using reversed vertical layers: from TOA to surface
    vertically_gridded_AK_da = xr.DataArray(vertically_gridded_AK,
                                            dims=("lat", "lon","MUSICA_layer"),
                                            coords={
                                                "lat": regionilats,
                                                "lon": regionilons,
                                                "MUSICA_layer": MUSICA_layers,
                                            })
    
    # return 
    return vertically_gridded_AK_da