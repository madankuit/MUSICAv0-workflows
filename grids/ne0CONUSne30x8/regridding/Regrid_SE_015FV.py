# This script includes functions used for CESM2.2 CAM-chem with SE dynamical core simulations
'''
Regrid_SE_015FV.py
this code is designed to regrid SE outputs (ne30 and MUSICA-V0) to regular lat and lon grids

Based on 'Rewrite_output.ipynb' from MUSICA Tutorial Nov. 2021

MODIFICATION HISTORY:
    Madankui Tao, 6, March, 2023: VERSION 1.00
    - Initial version
'''

### Module import ###
import numpy as np # for array manipulation and basic scientific calculation
import xarray as xr # To read NetCDF files
from scipy.interpolate import griddata # Simple regridding
from netCDF4 import Dataset # To write NetCDF files

#================================================================================================
def rewrite_ne30_1latlongrid_forCONUS_surflayer(diri,filename,varlist,savefilePath):
    """
    This function re-write the output from ne30 grid to a regular 1x1 degree lat and lon grid over CONUS
    Select lev=31 since the pressure system is positive down
    
        Input:
            Method: linear interpolation
            varlist: variables to re-grid
            savefilePath: path to the written .nc file (if == False would not save the output)
            
        Output:
            .nc file for the regrid dataset
        
        Use example:
            allvars_ds = rewrite_ne30_1latlongrid_forCONUS_surflayer(diri,filename,varlist,savefilePath)
        
    """
    # Define the lat/lon for desired grid 
    resolution = 1
    lon2d = np.arange(210,310+resolution,resolution)
    lat2d = np.arange(0,70+resolution,resolution)
    
    # This will put lat and lon into arrays
    X, Y = np.meshgrid(lon2d,lat2d)
    
    # empty dataset to store 
    allvars_ds = xr.Dataset()
    
    #-------------
    # get data on ncol, original output from the simulations
    ds_CONUS = xr.open_dataset(diri+filename)
    
    time_ar = ds_CONUS.time.values
    
    # Get the model lat/lon values for regridding
    mdllat = ds_CONUS['lat']
    mdllon = ds_CONUS['lon']
    
    for varidx in range(len(varlist)):
        vari = varlist[varidx]
        # regrid for this variable for the surface level; here the level using hybrid pressure positive down, so lev=31 is the surface layer
        vari_timearray = np.zeros([len(time_ar),len(lat2d),len(lon2d)])
        vari_timearray[:,:,:] = np.nan # initialize with nan
        
        # if there are multiple times
        for timeidx in range(len(time_ar)):
            vari_vals = griddata((mdllon,mdllat), ds_CONUS.isel(time=timeidx,lev=31)[vari], (X, Y), method='linear') #... using linear interpolation
            # add to the array
            vari_timearray[timeidx,:,:] = vari_vals
        
        # write to dataset for vari
        regrid_ds = xr.Dataset({
                    vari: xr.DataArray(
                                data   = vari_timearray,
                                dims   = ['time','lat','lon'],
                                coords = {'time':time_ar,'lat': lat2d,'lon': lon2d},
                                attrs  = {
                                    'units': ds_CONUS[vari].units, # stay the same in case vari is not concentration
                                    'long_name': ds_CONUS[vari].long_name
                                        }
                                    )
                                },
                        attrs = {'description': 'near surface estimated using linear interpolation'}
                    )
        # merge for all var
        if varidx==0:
            allvars_ds = regrid_ds
        else:
            # merge
            allvars_ds = xr.merge([allvars_ds, regrid_ds], compat='override')
            
    # save
    if savefilePath!=False:
        fullpathname_OUT = savefilePath+'CONUSRegrid1_'+filename #+'.nc'
        allvars_ds.to_netcdf(fullpathname_OUT)
        print("save to:",fullpathname_OUT)
        
    return allvars_ds

#================================================================================================
def rewrite_ne0CONUSne30x8_0125latlongrid_forCONUS_surflayer_v2(ds_CONUS,varlist,savefilePath,filename):
    """
    This function re-write the output from ne0CONUSne30x8 grid to a regular 0.125 x 0.125 lat and lon grid 
    Select lev=31 since the pressure system is positive down
    
        Input:
            Method: linear interpolation
            varlist: variables to re-grid
            savefilePath: path to the written .nc file (if == False would not save the output)
            filename: name of the saved file
            
        Output:
            .nc file for the regrid dataset
        
        Use example:
            ds_CONUS (on ne0CONUS30x8 grid)
            varlist = ['O3','CO','OH']
            savefilePath = '/home/taoma528/Scripts/CESM_analysis/Regrid_CONUS_data/'
            filename = 
            allvars_ds = rewrite_ne0CONUSne30x8_0125latlongrid_forCONUS_surflayer_v2(ds_CONUS,varlist,savefilePath)
        
    """
    # Define the lat/lon for desired grid 
    resolution = 0.125
    lat2d = np.arange(0,70+resolution,resolution)
    lon2d = np.arange(210,310+resolution,resolution)
    
    # This will put lat and lon into arrays
    X, Y = np.meshgrid(lon2d,lat2d)
    
    # empty dataset to store 
    allvars_ds = xr.Dataset()
    
    #-------------
    
    time_ar = ds_CONUS.time.values
    
    # Get the model lat/lon values for regridding
    mdllat = ds_CONUS['lat']
    mdllon = ds_CONUS['lon']
    
    for varidx in range(len(varlist)):
        vari = varlist[varidx]
        # regrid for this variable for the surface level; here the level using hybrid pressure positive down, so lev=31 is the surface layer
        vari_timearray = np.zeros([len(time_ar),len(lat2d),len(lon2d)])
        vari_timearray[:,:,:] = np.nan # initialize with nan
        
        # if there are multiple times
        for timeidx in range(len(time_ar)):
            vari_vals = griddata((mdllon,mdllat), ds_CONUS.isel(time=timeidx,lev=31)[vari], (X, Y), method='linear') #... using linear interpolation
            # add to the array
            vari_timearray[timeidx,:,:] = vari_vals
        
        # write to dataset for vari
        regrid_ds = xr.Dataset({
                    vari: xr.DataArray(
                                data   = vari_timearray,
                                dims   = ['time','lat','lon'],
                                coords = {'time':time_ar,'lat': lat2d,'lon': lon2d},
                                attrs  = {
                                    'units': ds_CONUS[vari].units,
                                    'long_name': ds_CONUS[vari].long_name
                                        }
                                    )
                                },
                        attrs = {'description': 'near surface estimated using linear interpolation'}
                    )
        # merge for all var
        if varidx==0:
            allvars_ds = regrid_ds
        else:
            # merge
            allvars_ds = xr.merge([allvars_ds, regrid_ds], compat='override')
            
    # save
    if savefilePath!=False:
        fullpathname_OUT = savefilePath+'CONUSRegrid0125_'+filename+'.nc'
        allvars_ds.to_netcdf(fullpathname_OUT)
        print("save to:",fullpathname_OUT)
        
    return allvars_ds

#================================================================================================
def rewrite_ne0CONUSne30x8_0125latlongrid_forCONUS_surflayer(diri,filename,varlist,savefilePath):
    """
    This function re-write the output from ne0CONUSne30x8 grid to a regular 0.125 x 0.125 lat and lon grid 
    Select lev=31 since the pressure system is positive down
    
        Input:
            Method: linear interpolation
            varlist: variables to re-grid
            savefilePath: path to the written .nc file (if == False would not save the output)
            
        Output:
            .nc file for the regrid dataset
        
        Use example:
            diri = EgFile_diri
            filename = 'f.e22.FCcotagsNudged.ne0CONUSne30x8.cesm220.2012-01.cam.h1.2013-08.nc'
            varlist = ['O3','CO','OH']
            savefilePath = '/home/taoma528/Scripts/CESM_analysis/Regrid_CONUS_data/'
            allvars_ds = rewrite_ne0CONUSne30x8_0125latlongrid_forCONUS_surflayer(diri,filename,varlist,savefilePath)
        
    """
    # Define the lat/lon for desired grid 
    resolution = 0.125
    lat2d = np.arange(0,70+resolution,resolution)
    lon2d = np.arange(210,310+resolution,resolution)
    
    # This will put lat and lon into arrays
    X, Y = np.meshgrid(lon2d,lat2d)
    
    # empty dataset to store 
    allvars_ds = xr.Dataset()
    
    #-------------
    # get data on ncol, original output from the simulations
    ds_CONUS = xr.open_dataset(diri+filename)
    
    time_ar = ds_CONUS.time.values
    
    # Get the model lat/lon values for regridding
    mdllat = ds_CONUS['lat']
    mdllon = ds_CONUS['lon']
    
    for varidx in range(len(varlist)):
        vari = varlist[varidx]
        # regrid for this variable for the surface level; here the level using hybrid pressure positive down, so lev=31 is the surface layer
        vari_timearray = np.zeros([len(time_ar),len(lat2d),len(lon2d)])
        vari_timearray[:,:,:] = np.nan # initialize with nan
        
        # if there are multiple times
        for timeidx in range(len(time_ar)):
            vari_vals = griddata((mdllon,mdllat), ds_CONUS.isel(time=timeidx,lev=31)[vari], (X, Y), method='linear') #... using linear interpolation
            # add to the array
            vari_timearray[timeidx,:,:] = vari_vals
        
        # write to dataset for vari
        regrid_ds = xr.Dataset({
                    vari: xr.DataArray(
                                data   = vari_timearray,
                                dims   = ['time','lat','lon'],
                                coords = {'time':time_ar,'lat': lat2d,'lon': lon2d},
                                attrs  = {
                                    'units': ds_CONUS[vari].units,
                                    'long_name': ds_CONUS[vari].long_name
                                        }
                                    )
                                },
                        attrs = {'description': 'near surface estimated using linear interpolation'}
                    )
        # merge for all var
        if varidx==0:
            allvars_ds = regrid_ds
        else:
            # merge
            allvars_ds = xr.merge([allvars_ds, regrid_ds], compat='override')
            
    # save
    if savefilePath!=False:
        fullpathname_OUT = savefilePath+'CONUSRegrid0125_'+filename+'.nc'
        allvars_ds.to_netcdf(fullpathname_OUT)
        print("save to:",fullpathname_OUT)
        
    return allvars_ds
