'''
This function is used to re-grid TROPOMI AK already processed to 0.15x0.15 degree to the vertical resolution of MUSICA for a given region on datei (YYYYMMDD)

# Revised from func_VerticalRegrid_TROPOMIAK_toMUSICAlevs.py, modified given different ways of processing MUSICA and TROPOMI

MODIFICATION HISTORY:
    13, December, 2023: VERSION 1.0
    - Initial version
    16, December, 2023: VERSION 1.1
    - Adjust for different casename
'''
# ============================================================
# CONFIGURATION - every path comes from config/paths.py, the single
# source of truth. Override cluster locations with the MUSICA_ENV_*
# environment variables documented there. Do not hard-code paths here.
# ============================================================
import sys as _sys, pathlib as _pathlib
_ROOT = next(_p for _p in _pathlib.Path(__file__).resolve().parents
             if (_p / 'config' / 'paths.py').exists())
_sys.path.insert(0, str(_ROOT))
from config.paths import (
    TROPOMI_015_DIRS, TROPOMI_L2_NO2_DIR, TROPOMI_L2_CO_DIR, TM5_SURFACE_P_DIR,
    REGRID_SUBDIR_BYVAR, case_hist_dir, case_regrid_dir,
)
# ============================================================


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

# Regridded TROPOMI L3 products, 0.15 deg (from config/paths.py)
L3_015TROPOMI_diri_dic = {k: str(v) + '/' for k, v in TROPOMI_015_DIRS.items()}

# Regridded TROPOMI air-mass-factor (AMF) fields, on the SAME 0.15° grid as the AK above.
# Needed for NO2 ONLY: the AK stored in the L2 PRODUCT group (and therefore in the regridded
# AK file read below) is the *total-column* averaging kernel. The tropospheric averaging
# kernel used in the tropospheric-VCD comparison is
#     AK_trop(l) = AK_total(l) × (AMF_total / AMF_trop)
# (TROPOMI ATBD S5P-KNMI-L2-0005-RP / Product User Manual Eq. 4). This mirrors the L2
# match/recalc convention (see tempo-v04-satellite-intercomparison, match_L2Scan.py
# _apply_ak_trop_conversion). HCHO is a (tropospheric) total-column retrieval, so its
# stored AK is used directly — no scaling. CO uses its total-column AK directly.
#
# Expected inputs, produced by the L2->L3 regrid step alongside the AK, each holding one
# 3-D variable 'TROPOMI_NO2_AMFtotal' / 'TROPOMI_NO2_AMFtrop' on (time, lat, lon):
#     <dir>S5P_RPRO_L2__NO2____AMFtotal_Global_Regrid015deg_<YYYYMMDD>.nc
#     <dir>S5P_RPRO_L2__NO2____AMFtrop_Global_Regrid015deg_<YYYYMMDD>.nc
L3_015AMF_diri_dic = {'NO2': str(TROPOMI_015_DIRS['NO2']) + '/'}


def _load_AMF_total_over_trop_ratio_015(varname, datei, lat_bot, lat_up, lon_right, lon_left):
    """NO2 only -- per-pixel AMF_total/AMF_trop ratio on the regridded 0.15deg grid.

    Used to convert the stored total-column averaging kernel to a tropospheric
    averaging kernel: AK_trop = AK_total x (AMF_total / AMF_trop). The ratio is
    layer-independent, so it is applied as a single scalar per (lat, lon) pixel.
    Pixels with AMF_trop <= 0 (or missing) are set to NaN. Returns a 2-D array
    aligned to the AK region grid (lat, lon).
    """
    amf_diri = L3_015AMF_diri_dic[varname]
    header   = fileheader_dic[varname]
    amf_total = xr.open_dataset(amf_diri+header+'AMFtotal_Global_Regrid015deg_'+datei+'.nc')['TROPOMI_'+varname+'_AMFtotal'].sel(time=datei)
    amf_trop  = xr.open_dataset(amf_diri+header+'AMFtrop_Global_Regrid015deg_'+datei+'.nc')['TROPOMI_'+varname+'_AMFtrop'].sel(time=datei)
    amf_total = amf_total.sel(lat=slice(lat_bot, lat_up), lon=slice(lon_right, lon_left)).values
    amf_trop  = amf_trop.sel(lat=slice(lat_bot, lat_up), lon=slice(lon_right, lon_left)).values
    with np.errstate(invalid='ignore', divide='ignore'):
        ratio = np.where(amf_trop > 0, amf_total / amf_trop, np.nan)
    return ratio

#=====Example use======
#--------------------------------------------------------------------------
# # Example use
# datei = '20180701'
# fileregion = 'Northeast'
# varname = 'HCHO' # 'HCHO';'NO2'

# vertically_gridded_AK_da = func_VerticalRegrid_TROPOMIAK_toMUSICAlevs_015latlon(datei, fileregion, varname)

#=====Main Function======
#--------------------------------------------------------------------------

def func_VerticalRegrid_TROPOMIAK_toMUSICAlevs_015latlon(casename, datei, fileregion, varname):
    """
    Calculate vertically gridded Average Kernel (AK) data to match MUSICA output.

    Parameters:
        datei (str): The date to process the Average Kernel.
        fileregion (str): The region of interest.
        varname (str): The variable name, either 'HCHO' or 'NO2'.

    Returns:
        xarray.DataArray: The vertically gridded Average Kernel data for the given date and region.
    """
    #--------------------------------------------------------------------------
    # Get the target grid
    lon_right,lat_bot,lon_left,lat_up = latlonbound(fileregion)
    
    #--------------------------------------------------------------------------
    # TROPOMI Test file for tm5_constant_a,tm5_constant_b,layers
    TROPOMI_L2NO2diri = str(TROPOMI_L2_NO2_DIR) + '/'
    testno2file = 'S5P_RPRO_L2__NO2____20180701T225858_20180702T004028_03711_03_020400_20221105T050222.nc'
    TROPOMI_NO2_PRODUCT = xr.open_dataset(TROPOMI_L2NO2diri+testno2file,group="/PRODUCT/",engine="netcdf4")
    TROP_NO2_egda = xr.open_dataset(TROPOMI_L2NO2diri+testno2file,group="PRODUCT/SUPPORT_DATA/INPUT_DATA/").sel(time=0)
    # TM5 constant for a and b are the same for all files
    tm5_constant_a = TROPOMI_NO2_PRODUCT.tm5_constant_a
    tm5_constant_b = TROPOMI_NO2_PRODUCT.tm5_constant_b
    layers = TROPOMI_NO2_PRODUCT.layer.values
    layernum = len(layers)    
    
    #--------------------------------------------------------------------------
    # Read in test MUSICA output file: get the vertical levels from MUSICA; does not matter which specific case
    # Representative MUSICA h2 file, used only to read the vertical levels.
    test_h2_MUSICAfile = str(case_hist_dir(casename) / f'{casename}.cam.h2.2018-07-01-03600.nc')
    testh2_ds = xr.open_dataset(test_h2_MUSICAfile)
    # hybrid level at midpoints
    MUSICA_layers = testh2_ds.lev.values
    MUSICAlaynum = len(MUSICA_layers)
    del test_h2_MUSICAfile,testh2_ds

    # PMID processed to 0.15 degree already for the given case
    regridByVar_diri = str(case_regrid_dir(casename, REGRID_SUBDIR_BYVAR)) + '/'
    MUSICA_PMID = xr.open_dataset( regridByVar_diri+casename+'.cam.h2.MassConserve_latlon015.PMID.20180701T20180801.nc' )
    regioni_MUSICA_PMID = MUSICA_PMID.PMID.sel(lat=slice(lat_bot, lat_up),lon=slice(lon_right,lon_left))
    # Approx at TROPOMI overpass time
    dateFMT2 = datei[:4]+'-'+datei[4:6]+'-'+datei[6:8]
    averageTimeStart,averageTimeEnd = timezone_approx1330LT(fileregion)
    sel_times = [np.datetime64(dateFMT2+averageTimeStart),np.datetime64(dateFMT2+averageTimeEnd)]
    # select for the approx two hours at 1 and 2 PM, then take the average
    regioni_MUSICA_PMID_1330LT = regioni_MUSICA_PMID.sel(time=sel_times).mean(dim='time') 
    
    #--------------------------------------------------------------------------
    # Calculate TM5 PMID for the given fileregion and datei, already processed to 0.15 degree lat-lon
    # TM5PMID_datei_layers have layers revserved from TOA to surface
    TM5SurfaceP_FileOut_diri = str(TM5_SURFACE_P_DIR) + '/'
    TM5SurfaceP_Path = TM5SurfaceP_FileOut_diri+'TM5MP_Model_SurfaceP_Global_015deg_2018-06-30T2018-08-02.nc'
    TM5SurfaceP_ds = xr.open_dataset(TM5SurfaceP_Path)
    regioniSurfaceP_ds = TM5SurfaceP_ds.sel(lat=slice(lat_bot,lat_up),lon=slice(lon_right,lon_left))
    dateiPs = regioniSurfaceP_ds.Ps.sel(time=datei).values
    # get the horizontal resolution grid
    regionilats = regioniSurfaceP_ds.lat
    regionilons = regioniSurfaceP_ds.lon

    # Calculate the absolute pressure in hPa | use the pressure at the middle of the interface level (match MUSICA) 
    TM5PMID_datei_layers = np.zeros([len(regionilats),len(regionilons),layernum])
    TM5PMID_datei_layers[:,:,:] = np.nan
    # calculate the pressure levels PMID using Ps_alldays_data from the last layer to use for the reversed 
    for LAYidx, reversedLAYidx in enumerate(range(layernum - 1, -1, -1)):
        # print(LAYidx,reversedLAYidx)
        # pressure at the lowest point [0] for lower and [1] for upper interface level
        tm5_constant_a_layi = np.nanmean(tm5_constant_a[reversedLAYidx].values)
        tm5_constant_b_layi = np.nanmean(tm5_constant_b[reversedLAYidx].values)

        TM5PMID_layi = tm5_constant_a_layi+tm5_constant_b_layi*(dateiPs)
        # put the revsered layer to the layer count from beginning
        TM5PMID_datei_layers[:,:,LAYidx] = TM5PMID_layi
        
    # Save using reversed vertical layers
    TM5PMID_datei_ds = xr.DataArray(TM5PMID_datei_layers,
                                            dims=("lat", "lon","reversed_layer"),
                                            coords={
                                                "lat": regionilats,
                                                "lon": regionilons,
                                                "reversed_layer": layers,
                                            })
    
    #--------------------------------------------------------------------------
    # Read in TROPOMI AK for the given fileregion and datei; reverse the vertical layer order
    AK_FileOutPath = L3_015TROPOMI_diri_dic[varname]+fileheader_dic[varname]+'AK_Global_Regrid015deg_'+datei+'.nc'
    AK_ds = xr.open_dataset(AK_FileOutPath)
    # select for the date
    AK_xar = AK_ds['TROPOMI_'+varname+'_AK'].sel(time=datei)

    # Reverse the 'layer' dimension indices in TROPOMI AK
    AK_xar_reversed = AK_xar.isel(layer=slice(None, None, -1))
    # Reverse the data values along the 'layer' dimension
    AK_xar_reversed['layer'] = AK_xar_reversed['layer'].values[::-1]
    # You can also rename the 'layer' dimension to 'reversed_layer' if needed
    AK_xar_reversed = AK_xar_reversed.rename({'layer': 'reversed_layer'})
    # select for the given region
    regioni_datei_AK_xar_reversed = AK_xar_reversed.sel(lat=slice(lat_bot, lat_up),lon=slice(lon_right,lon_left))
    
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
            values_source = regioni_datei_AK_xar_reversed.sel(lat=lati,lon=loni).values
            # Pressure levels on the source grid, used by TROPOMI in Pa
            pressure_source = TM5PMID_datei_ds.sel(lat=lati,lon=loni).values
            # Pressure levels on the target grid, used PMID by MUSICA in Pa
            pressure_target = regioni_MUSICA_PMID_1330LT.sel(lat=lati,lon=loni).values
            # Use numpy.interp to interpolate values from the source grid to the target grid
            values_target = np.interp(pressure_target, pressure_source, values_source)
            vertically_gridded_AK[latidx,lonidx,:] = values_target
            # print("Interpolated values on the target grid:", values_target)

    #--------------------------------------------------------------------------
    # NO2 only: convert the stored total-column AK to a tropospheric AK using the
    # per-pixel AMF ratio, AK_trop = AK_total x (AMF_total / AMF_trop). The ratio is
    # vertically uniform, so multiply every MUSICA layer by the 2-D (lat, lon) field.
    # HCHO (tropospheric total-column retrieval) uses its stored AK directly.
    if varname == 'NO2':
        amf_ratio = _load_AMF_total_over_trop_ratio_015(varname, datei, lat_bot, lat_up, lon_right, lon_left)
        vertically_gridded_AK = vertically_gridded_AK * amf_ratio[:, :, np.newaxis]

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