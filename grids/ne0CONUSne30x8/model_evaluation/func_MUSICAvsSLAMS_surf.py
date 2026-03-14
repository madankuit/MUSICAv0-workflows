'''
MUSICAvsSLAMS_surf.py
the functions are designed to compare MUSICA model simulations with SLAMS/AQS observations
MODIFICATION HISTORY:
    Madankui Tao, 3, OCT, 2023: VERSION 1.0
    - Initial version
    Madankui Tao, 18, DEC, 2023: VERSION 1.1
    - Add a modified function for hourly AQS files
'''

#================================================================================================
# ### functions import ###
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
import sys
sys.path.append('/home/taoma528/Scripts/CESM_analysis/functions')
from func_ModelEval_statistical_tests import *
from SE_analysis import get_site_index

# Grid
EgFile_diri = '/home/taoma528/Scripts/CESM_analysis/tutorial_files/'
# Read SCRIP file that has grid information needed to plot values on a map
SCRIP_CONUS = EgFile_diri+'ne0CONUS_ne30x8_np4_SCRIP.nc'
SCRIP_ne30 = EgFile_diri+'ne30np4_091226_pentagons.nc'

### preset variables ###
lev_idx = -1 # -1 for vertical level 

#================================================================================================
# Version 0 function
### Approx column index based on (lat,lon) ###
def write_MonitorIDtoMUSICAcolidx_csvfile(varname,AQSvar_df):
    """
    The function calculates the approximate column index of the MUSICA grid for each unique MonitorID location and saves the results to a CSV file.
    Inputs:
        varname (str): The name of the species or variable for which column indices are being calculated.
        AQSvar_df (pandas.DataFrame): A DataFrame containing data for the specified species, including 'MonitorID', 'Latitude', and 'Longitude'.

    CSVs are saved in the 'colidxCSV' directory with filenames like 'MUSICA0_AQS{varname}monitor_colidx.csv'.
    Each CSV contains columns:
    - 'MonitorID': Unique MonitorID values
    - 'Latitude': Latitude values for MonitorID locations
    - 'Longitude': Longitude values for MonitorID locations
    - 'MUSICA0_colIndex': Approximate column index on the MUSICA grid
    - 'Approx_MUSICA0_lat': Approximate latitude on the MUSICA grid
    - 'Approx_MUSICA0_lon': Approximate longitude on the MUSICA grid
    """
    
    # location for the output
    CVSOUT_diri = '/home/taoma528/Scripts/CESM_analysis/colidxCSV/'
    
    # Read an example MUSICA ne0CONUSne30x8 ds for model lat and lon
    SpinUp_diri = '/net/fs09/d0/taoma528/CESM22/archive/f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.base.Y2018SpinupCmip6hist2013/atm/hist/'
    ex_h1_filename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.base.Y2018SpinupCmip6hist2013.cam.h1.2018-05-01-00000.nc'
    ds_CONUS = xr.open_dataset(SpinUp_diri+ex_h1_filename)
    
    # Keep only the unique MonitorID values along with their Latitude and Longitude
    unique_AQSmonitor_locations = AQSvar_df[['MonitorID', 'Latitude', 'Longitude']].drop_duplicates().reset_index(drop=True)
    
    # loop through the monitors and find the corresponding column index
    ls_Index_MonitorIDi = []
    ls_Model_lati = []
    ls_Model_loni = []

    for MonitorIDi in unique_AQSmonitor_locations.MonitorID.values:
        # Get the latitude and longitude values for 'MonitorIDi'
        MonitorIDi_df = unique_AQSmonitor_locations[unique_AQSmonitor_locations['MonitorID'] == MonitorIDi]

        lati = MonitorIDi_df['Latitude'].iloc[0]
        loni = MonitorIDi_df['Longitude'].iloc[0]

        # find the model index
        Index_MonitorIDi = get_site_index( site_lat=lati, site_lon=360+loni, scrip_file=SCRIP_CONUS )
        if Index_MonitorIDi==None:
            # add to the list
            ls_Index_MonitorIDi.append('Find None')
            ls_Model_lati.append(Model_lati)
            ls_Model_loni.append(Model_loni)

        else:
            Model_lati = ds_CONUS.lat.values[Index_MonitorIDi]
            Model_loni = ds_CONUS.lon.values[Index_MonitorIDi]
            # add to the list
            ls_Index_MonitorIDi.append(Index_MonitorIDi)
            ls_Model_lati.append(Model_lati)
            ls_Model_loni.append(Model_loni)

    # append to the df
    unique_AQSmonitor_locations['MUSICA0_colIndex'] = ls_Index_MonitorIDi
    unique_AQSmonitor_locations['Approx_MUSICA0_lat'] = ls_Model_lati
    unique_AQSmonitor_locations['Approx_MUSICA0_lon'] = ls_Model_loni

    # Save the DataFrame to a CSV file
    file_path = CVSOUT_diri+'MUSICA0_AQS'+varname+'Monitor_colidx.csv'
    unique_AQSmonitor_locations.to_csv(file_path, index=False)
    
    print("Saved to:",file_path)
    
#================================================================================================
# Modified for hourly dataframe
def write_hourlyAQS_MonitorIDtoMUSICAcolidx_csvfile(varname,AQSvar_df,CVSOUT_diri):
    """
    The function calculates the approximate column index of the MUSICA grid for each unique MonitorID location and saves the results to a CSV file.
    Inputs:
        varname (str): The name of the species or variable for which column indices are being calculated.
        AQSvar_df (pandas.DataFrame): A DataFrame containing data for the specified species, including 'MonitorID', 'Latitude', and 'Longitude'.
        CVSOUT_diri: a path to a directory to save

    CSVs are saved in the 'colidxCSV' directory with filenames like MUSICA_index_diri+'MUSICA0_hourlyAQS_{varname}_Monitor_colidx.csv'
    Each CSV contains columns:
    - 'MonitorID': Unique MonitorID values
    - 'Latitude': Latitude values for MonitorID locations
    - 'Longitude': Longitude values for MonitorID locations
    - 'MUSICA0_colIndex': Approximate column index on the MUSICA grid
    - 'Approx_MUSICA0_lat': Approximate latitude on the MUSICA grid
    - 'Approx_MUSICA0_lon': Approximate longitude on the MUSICA grid
    """
    
    # Read an example MUSICA ne0CONUSne30x8 ds for model lat and lon
    SpinUp_diri = '/net/fs09/d0/taoma528/CESM22/archive/f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.base.Y2018SpinupCmip6hist2013/atm/hist/'
    ex_h1_filename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.base.Y2018SpinupCmip6hist2013.cam.h1.2018-05-01-00000.nc'
    ds_CONUS = xr.open_dataset(SpinUp_diri+ex_h1_filename)
    
    # Keep only the unique MonitorID values along with their Latitude and Longitude
    unique_AQSmonitor_locations = AQSvar_df[['MonitorID', 'Latitude', 'Longitude']].drop_duplicates().reset_index(drop=True)
    
    # loop through the monitors and find the corresponding column index
    ls_Index_MonitorIDi = []
    ls_Model_lati = []
    ls_Model_loni = []

    for MonitorIDi in unique_AQSmonitor_locations.MonitorID.values:
        # Get the latitude and longitude values for 'MonitorIDi'
        MonitorIDi_df = unique_AQSmonitor_locations[unique_AQSmonitor_locations['MonitorID'] == MonitorIDi]

        lati = MonitorIDi_df['Latitude'].iloc[0]
        loni = MonitorIDi_df['Longitude'].iloc[0]

        # find the model index
        Index_MonitorIDi = get_site_index( site_lat=lati, site_lon=360+loni, scrip_file=SCRIP_CONUS )
        if Index_MonitorIDi==None:
            # add to the list
            ls_Index_MonitorIDi.append('Find None')
            ls_Model_lati.append(Model_lati)
            ls_Model_loni.append(Model_loni)

        else:
            Model_lati = ds_CONUS.lat.values[Index_MonitorIDi]
            Model_loni = ds_CONUS.lon.values[Index_MonitorIDi]
            # add to the list
            ls_Index_MonitorIDi.append(Index_MonitorIDi)
            ls_Model_lati.append(Model_lati)
            ls_Model_loni.append(Model_loni)

    # append to the df
    unique_AQSmonitor_locations['MUSICA0_colIndex'] = ls_Index_MonitorIDi
    unique_AQSmonitor_locations['Approx_MUSICA0_lat'] = ls_Model_lati
    unique_AQSmonitor_locations['Approx_MUSICA0_lon'] = ls_Model_loni

    # Save the DataFrame to a CSV file
    file_path = CVSOUT_diri+'MUSICA0_hourlyAQS_'+varname+'_Monitor_colidx.csv'
    unique_AQSmonitor_locations.to_csv(file_path, index=False)
    
    print("Saved to:",file_path)


#================================================================================================
# ### h1 -daily file | Scatter Plot of Model vs. Observation ###
### scatter plot with statistical test ###
def plot_scatter_with_regression(SurfObs_ar, MUSICA_ar, Title, varname, axislim, givenpval, numbertype):
    """
    Create a scatter plot with regression line.

    Parameters:
        SurfObs_ar (array-like): Array of observed values.
        MUSICA_ar (array-like): Array of MUSICA model values.
        Title (str): Title for the plot.
        varname (str): name of the chemical species
        yaxislim (float or int): max lim of x-axis and y-aixs
        
        # the names are what's used in MUSICA modeling

    Returns:
        None
    """
    varlabel_dic = {'O3':r'$O_{3}$','NO2':r'$NO_{2}$',
                    'CO':r'$CO$','SO2':r'$SO_{2}$',
                    'PM25':r'$PM_{2.5}$',
                    'T':r'T',
                    'C2H2':r'$C_{2}H_{2}$','C2H6':r'$C_{2}H_{6}$',
                    'CH2O':r'$CH_{2}O$','C3H6':r'$C_{3}H_{6}$',
                    'C3H8':r'$C_{3}H_{8}$','ISOP':r'$Isoprene$',
                    'BENZENE':r'$Benzene$'}
    varunit_dic = {'O3':'ppb','NO2':'ppb',
                   'CO':'ppb','SO2':'ppb',
                   'PM25':r'$µg/m^{3}$',
                    'T':'Temperature',
                   }
    
    # Calculate the R2,rma_slope,and rma_intercept of the regression line.
    R2,rma_slope,rma_intercept = model_vs_obs_RMA(SurfObs_ar,MUSICA_ar, givenpval, numbertype)
    
    fig, ax = plt.subplots(figsize=(8, 6))
    plt.grid(b=None)

    # Text to display R-squared and slope
    textstr = (r'$r^{2}$=%.2f; slope=%.2f' % (R2, rma_slope))

    # Create the scatter plot
    sc = ax.scatter(SurfObs_ar, MUSICA_ar, c='black',
                    vmin=0, vmax=80, s=20, marker='o', label=textstr)

    # Set labels and title
    ax.set_xlabel(r'SLAMS '+varlabel_dic[varname]+' ['+varunit_dic[varname]+']', fontsize=22)
    ax.set_ylabel(r'MUSICA '+varlabel_dic[varname]+' ['+varunit_dic[varname]+']', fontsize=22)
    ax.set_title(Title, y=1.03, fontsize=22)
    
    # Set tick label font size
    ax.tick_params(axis='both', which='major', labelsize=16)

    # Set plot limits
    ax.set_xlim([0, axislim])
    ax.set_ylim([0, axislim])

    # Plot 1:1 line in red
    ax.plot(SurfObs_ar, SurfObs_ar, 'r-', label="1:1 line")

    # Calculate and plot the regression line
    predicted_values = rma_slope * SurfObs_ar + rma_intercept
    plt.plot(SurfObs_ar, predicted_values, color='blue', label='Line of Best-fit')

    # Add legend
    plt.legend(loc="lower right", fontsize=16, markerscale=1)

    # Show the plot
    plt.show()

# Example usage:
# plot_scatter_with_regression(SurfObs_ar, MUSICA_ar, datei)

#================================================================================================
# ### h1 -daily file | Map of Model minus Observation ###
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

def plot_difference_on_CONUSmap(latobs, lonobs, valobs, rangemax, Title):
    """
    Create a scatter plot on a map of the Continental United States (CONUS) to visualize differences.

    Parameters:
        latobs (array-like): Array of latitudes.
        lonobs (array-like): Array of longitudes.
        valobs (array-like): Array of values to visualize.
        rangemax (float): Maximum range for the color scale.
        Title (str): Title for the plot.

    Returns:
        None
    """
    # Create a map using PlateCarree projection
    fig, ax = plt.subplots(figsize=(14, 10), subplot_kw={'projection': ccrs.PlateCarree()})

    # Set the extent to the desired lon_range and lat_range
    lon_range = [-140, -50]
    lat_range = [15, 60]
    ax.set_extent([lon_range[0], lon_range[1], lat_range[0], lat_range[1]], crs=ccrs.PlateCarree())

    # Add state, country, and coastline boundaries as features
    ax.add_feature(cfeature.BORDERS, linestyle='-', linewidth=0.5)
    ax.add_feature(cfeature.STATES, linestyle='-', linewidth=0.5)
    ax.add_feature(cfeature.COASTLINE, linestyle='-', linewidth=0.5)

    # Plot the data as scattered points colored by value
    sc = ax.scatter(lonobs, latobs, c=valobs, cmap='RdBu_r', s=20,
                    vmin=-rangemax, vmax=rangemax, transform=ccrs.PlateCarree())

    # Add colorbar with specified fontsize, size, and position
    cbar = plt.colorbar(sc, label='[ppb]', ax=ax, location='bottom', pad=0.05, aspect=40, shrink=0.8)
    cbar.ax.tick_params(labelsize=14)  # Increase fontsize of tick labels
    cbar.ax.set_xlabel('[ppb]', fontsize=16)  # Set fontsize of colorbar label

    # Set fontsize for various labels
    ax.set_title(Title, y=1.02, fontsize=16)
    ax.set_xlabel('Longitude', fontsize=14)
    ax.set_ylabel('Latitude', fontsize=14)
    ax.tick_params(axis='both', labelsize=12)

    # Show the plot
    plt.show()

# Example usage:
# plot_difference_on_CONUSmap(latobs, lonobs, valobs, rangemax, Title)
#================================================================================================
# 



