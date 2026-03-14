# This function is used to re-grid TROPOMI AK already processed to 0.1x0.1 degree to the vertical resolution of MUSICA for a given region on datei (YYYYMMDD)
# Created by Madankui Tao on Nov 13, 2023

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

def latlonbound(fileregion):
    """This function returns the lats and lons boundary of a given region
        Regions include: West, Midwest, Southwest, Southeast, Northeast
    """
    # Define regions        
    if fileregion == "WestCoast":
        lon_right,lat_bot,lon_left,lat_up = [-125, 31, -115, 49]
    elif fileregion == "Mountain":
        lon_right,lat_bot,lon_left,lat_up = [-115, 31, -103, 49]
    elif fileregion == "Midwest":
        lon_right,lat_bot,lon_left,lat_up = [-103, 36.5, -85, 49]
    elif fileregion == "Southwest":
        lon_right,lat_bot,lon_left,lat_up = [-103, 25.5, -85, 36.5]
    elif fileregion == "Southeast":
        lon_right,lat_bot,lon_left,lat_up = [-85, 24.5, -75, 36.5]
    elif fileregion == "Northeast":
        lon_right,lat_bot,lon_left,lat_up = [-85, 36.5, -67, 49]
    
    return [lon_right,lat_bot,lon_left,lat_up ]

#=====Main Function======
#--------------------------------------------------------------------------
# # Example use
# datei = '20180702'

# # get TM5PMID_ds.TM5PMID for one day
# SurfaceP_FileOut_diri = "/net/fs09/d0/taoma528/Datasets/TROPOMI_RegridL2_ToL3/Global/TROPOMI_SurfaceP/"
# startDate = "2018-07-02"
# endDate = "2018-07-03"
# SurfaceP_Path = SurfaceP_FileOut_diri+'TROPOMI_SurfaceP_Global_01deg_'+startDate+'T'+endDate+'.nc'
# SurfaceP_ds = xr.open_dataset(SurfaceP_Path)

# fileregion = 'Northeast'

# varname = 'NO2' # HCHO

#=====Main Function======
#--------------------------------------------------------------------------

def func_VerticalRegrid_TROPOMIAK_toMUSICAlevs(datei, TM5SurfaceP_ds, fileregion, varname):
    """
    Calculate vertically gridded Average Kernel (AK) data to match MUSICA output.

    Parameters:
        datei (str): The date to process the Average Kernel.
        TM5SurfaceP_ds (xarray.Dataset): The surface pressure dataset from TROPOMI (using TM5 Model).
        fileregion (str): The region of interest.
        varname (str): The variable name, either 'HCHO' or 'NO2'.

    Returns:
        xarray.DataArray: The vertically gridded Average Kernel data for the given date and region.
    """
    #--------------------------------------------------------------------------
    # Get the target grid
    lon_right,lat_bot,lon_left,lat_up = latlonbound(fileregion)
    # One file as an example to get the 0.01x0.01 grid
    TROPOMI_var = xr.open_dataset('/net/fs09/d0/taoma528/Datasets/TROPOMI_RegridL2_ToL3/Global/TROPOMI_NO2_01deg/TROPOMI_NO2_Global_01deg_20180801.nc')
    # select for the domain
    region_TROPOMI_var = TROPOMI_var.sel(lat=slice(lat_bot, lat_up),lon=slice(lon_right,lon_left))
    regionilats01 = region_TROPOMI_var.lat.values
    regionilons01 = region_TROPOMI_var.lon.values
    
    #--------------------------------------------------------------------------
    # TROPOMI Test file for tm5_constant_a,tm5_constant_b,layers
    TROPOMI_L2NO2diri = "/net/fs09/d0/taoma528/Datasets/TROPOMI_NO2_L2_RPROv10202_2018_Global/"
    testno2file = 'S5P_RPRO_L2__NO2____20181015T060728_20181015T075055_05205_01_010202_20190227T185503.nc'
    TROPOMI_NO2_PRODUCT = xr.open_dataset(TROPOMI_L2NO2diri+testno2file,group="/PRODUCT/",engine="netcdf4")
    TROP_NO2_egda = xr.open_dataset(TROPOMI_L2NO2diri+testno2file,group="PRODUCT/SUPPORT_DATA/INPUT_DATA/").sel(time=0)

    Origin_latitude = TROPOMI_NO2_PRODUCT.latitude[0,:,:] #.values
    Origin_longitude = TROPOMI_NO2_PRODUCT.longitude[0,:,:]#.values
    tm5_constant_a = TROPOMI_NO2_PRODUCT.tm5_constant_a
    tm5_constant_b = TROPOMI_NO2_PRODUCT.tm5_constant_b
    layers = TROPOMI_NO2_PRODUCT.layer.values
    layernum = len(layers)
    
    #--------------------------------------------------------------------------
    # Test RegridL1 MUSICA file: get the vertical levels from MUSICA
    RegridL1_diri = '/net/fs09/d0/taoma528/CESM22/Regridded_MUSICA_Output/2018_1330LT_TROPOMIcomp/latlon01_MUSICAoutput/'
    ProcessLevel = 'L1.'
    # search for the file
    fileregion_filelist = glob.glob(os.path.join(RegridL1_diri, f'*{fileregion}*'))
    combined_file_list = [file for file in fileregion_filelist if all(keyword in file for keyword in [ProcessLevel])]
    if len(combined_file_list)==0:
        print('No MUSICA file on %s over %s'%(fileregion))
    else:
        dateifile = sorted(combined_file_list)[0]
    RegridL1_ds = xr.open_dataset(dateifile)
    regioni_MUSICA_PMID = RegridL1_ds.PMID.sel(lat=slice(lat_bot, lat_up),lon=slice(lon_right,lon_left))
    MUSICAlaynum = len(regioni_MUSICA_PMID.lev)
    MUSICA_layers = regioni_MUSICA_PMID.lev.values
    
    #--------------------------------------------------------------------------
    # Calculate TM5 PMID
    # TM5PMID_datei_layers have layers revserved from TOA to surface
    regioniSurfaceP_ds = TM5SurfaceP_ds.sel(lat=slice(lat_bot,lat_up),lon=slice(lon_right,lon_left))
    dateiPs = regioniSurfaceP_ds.Ps.sel(time=datei).values

    # Calculate the absolute pressure in hPa | use the pressure at the middle of the interface level (match MUSICA) 
    TM5PMID_datei_layers = np.zeros([len(regionilats01),len(regionilons01),layernum])
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
                                                "lat": regionilats01,
                                                "lon": regionilons01,
                                                "reversed_layer": layers,
                                            })
    #--------------------------------------------------------------------------
    # Read in TROPOMI AK
    if varname=='NO2':
        AK_FileOut_diri =  "/net/fs09/d0/taoma528/Datasets/TROPOMI_RegridL2_ToL3/Global/TROPOMI_NO2AK_01deg/"
        AK_FileOutPath = AK_FileOut_diri+'TROPOMI_NO2AK_Global_01deg_'+datei+'.nc'
        AK_ds = xr.open_dataset(AK_FileOutPath)
        AK_xar = AK_ds.TROPOMI_NO2_AK.sel(time=datei)
        
    elif varname=='HCHO':
        AK_FileOut_diri =  "/net/fs09/d0/taoma528/Datasets/TROPOMI_RegridL2_ToL3/Global/TROPOMI_HCHOAK_01deg/"
        AK_FileOutPath = AK_FileOut_diri+'TROPOMI_HCHOAK_Global_01deg_'+datei+'.nc'
        AK_ds = xr.open_dataset(AK_FileOutPath)
        AK_xar = AK_ds.TROPOMI_HCHO_AK.sel(time=datei)
        
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
    vertically_gridded_AK = np.zeros([len(regionilats01),len(regionilons01),MUSICAlaynum])

    # reverse the order for the vertical layers
    for latidx in range(len(regionilats01)):
        lati = regionilats01[latidx]
        for lonidx in range(len(regionilons01)):
            loni = regionilons01[lonidx]
            # print(latidx,lonidx)
            # Values on the source grid
            values_source = regioni_datei_AK_xar_reversed.sel(lat=lati,lon=loni).values
            # Pressure levels on the source grid, used by TROPOMI in Pa
            pressure_source = TM5PMID_datei_ds.sel(lat=lati,lon=loni).values
            # Pressure levels on the target grid, used PMID by MUSICA in Pa
            pressure_target = regioni_MUSICA_PMID.sel(lat=lati,lon=loni).values
            # Use numpy.interp to interpolate values from the source grid to the target grid
            values_target = np.interp(pressure_target, pressure_source, values_source)
            vertically_gridded_AK[latidx,lonidx,:] = values_target
            # print("Interpolated values on the target grid:", values_target)
    
    #--------------------------------------------------------------------------
    # Write to a xarray
    # Save using reversed vertical layers
    vertically_gridded_AK_da = xr.DataArray(vertically_gridded_AK,
                                            dims=("lat", "lon","MUSICA_layer"),
                                            coords={
                                                "lat": regionilats01,
                                                "lon": regionilons01,
                                                "MUSICA_layer": MUSICA_layers,
                                            })
    
    # return 
    return vertically_gridded_AK_da