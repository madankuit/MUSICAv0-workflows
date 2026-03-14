'''
ReadMUSICAOutput.py
this code is designed to read in MUSICA model simulation outputs 

MODIFICATION HISTORY:
    Madankui Tao, 3, OCT, 2023: VERSION 1.00
    - Initial version
    Madankui Tao, 19, FEB, 2024: VERSION 2.00
    - Incorporate the old 'July5Day' ipynb functions here
    
'''

### Module import ###
import numpy as np # for array manipulation and basic scientific calculation
import xarray as xr # To read NetCDF files
import matplotlib.pyplot as plt # Core library for plotting
import matplotlib.cm as cm # To use different colormaps

import warnings
# Ignore a specific warning
warnings.filterwarnings('ignore', message='Some warning message')
# Ignore the FutureWarning caused by iteritems()
warnings.filterwarnings("ignore", category=FutureWarning, message=".*iteritems.*")

import os
import glob

###########
# Set the path to your SCRIP grid files, e.g.:
# SCRIP_CONUS = '/path/to/ne0CONUS_ne30x8_np4_SCRIP.nc'
# SCRIP_ne30  = '/path/to/ne30np4_091226_pentagons.nc'

###########
import os
import glob
import xarray as xr

###########
def BasicRead_ds(RunPath,histfreq):
    """
    Read a daily dataset from a specific directory

    Parameters:
    RunPath (str): Path to the directory containing the files.
    
    Returns:
    xarray.DataArray: .
    """
    # Get a list of file paths in the specified directories with the given keyword
    Run_filelist = glob.glob(os.path.join(RunPath, f'*{histfreq}*'))

    # Pick the first file from the list (assuming only one file per directory)
    RunFile = Run_filelist[0]
    
    # Read in the dataset from the file
    Run_ds = xr.open_dataset(RunFile)
    print('Read in files %s' % (RunFile))
    
    return Run_ds

###########
def Read_histh1_Runds(RunPath, datei, lev_idx, varname):
    """
    Read a daily dataset from a specific directory and extract the variable data for a given time and vertical level.

    Parameters:
    RunPath (str): Path to the directory containing the files.
    time_idx (int): Index of the desired time slice.
    lev_idx (int): Index of the desired vertical level (e.g., pressure level).
    varname (str): Name of the variable to extract. if None, read the entire dataset

    Returns:
    xarray.DataArray: The variable data for the specified time and vertical level.
    """
    # for daily average file
    histfreq = 'h1'
    datei_str = datei+'T00:00:00.000000000'
    
    # Get a list of file paths in the specified directories with the given keyword
    Run_filelist = glob.glob(os.path.join(RunPath, f'*{histfreq}*'))

    # Pick the first file from the list (assuming only one file per directory)
    RunFile = Run_filelist[0]
    
    # Read in the dataset from the file
    Run_ds = xr.open_dataset(RunFile)
    print('Read in files %s' % (RunFile))
    
    # Select the variable data for the specified time and vertical level from the dataset
    # Run_DayLevi_Vari = Run_ds.isel(lev=lev_idx)
    Run_DayLevi_Vari = Run_ds.isel(lev=lev_idx).sel(time=datei_str)
    
    if varname!='None':
        Run_DayLevi_Vari = Run_DayLevi_Vari[varname]
        
    print('Run - Time: %s; Vertical Level: %s' % (str(Run_DayLevi_Vari.time.values),
                                                 str(int(Run_DayLevi_Vari.lev.values)) + ' hPa'))
    
    return Run_DayLevi_Vari

########### 
def Read_histh1_Runds_v2(RunPath, date_idx, lev_idx, varname):
    """
    Read a daily dataset from a specific directory and extract the variable data for a given time and vertical level.

    Parameters:
    RunPath (str): Path to the directory containing the files.
    time_idx (int): Index of the desired time slice.
    lev_idx (int): Index of the desired vertical level (e.g., pressure level).
    varname (str): Name of the variable to extract.

    Returns:
    xarray.DataArray: The variable data for the specified time and vertical level.
    """
    # for daily average file
    histfreq = 'h1'
    
    # Get a list of file paths in the specified directories with the given keyword
    Run_filelist = glob.glob(os.path.join(RunPath, f'*{histfreq}*'))

    # Pick the first file from the list (assuming only one file per directory)
    RunFile = Run_filelist[0]
    
    # Read in the dataset from the file
    Run_ds = xr.open_dataset(RunFile)
    print('Read in files %s' % (RunFile))
    
    # Select the variable data for the specified time and vertical level from the dataset
    # Run_DayLevi_Vari = Run_ds.isel(lev=lev_idx)
    Run_DayLevi_Vari = Run_ds.isel(lev=lev_idx,time=date_idx)
    
    if varname!='None':
        Run_DayLevi_Vari = Run_DayLevi_Vari[varname]
        
    print('Run - Time: %s; Vertical Level: %s' % (str(Run_DayLevi_Vari.time.values),
                                                 str(int(Run_DayLevi_Vari.lev.values)) + ' hPa'))
    
    return Run_DayLevi_Vari
########### 
def Read_histh2_Runds(RunPath, datei, lev_idx, varname):
    """
    Read a hourly dataset from a specific directory and extract the variable data for a given date and vertical level.

    Parameters:
    RunPath (str): Path to the directory containing the files.
    datei (str): The desired time slice.
    lev_idx (int): Index of the desired vertical level (e.g., pressure level).
    varname (str): Name of the variable to extract.

    Returns:
    xarray.DataArray: The variable data for the specified time and vertical level.
    """
    # for hourly average file
    histfreq = 'h2'
    searchstr = histfreq + '.' + datei
    
    # Get a list of file paths in the specified directories with the given keyword
    Run_filelist = glob.glob(os.path.join(RunPath, f'*{searchstr}*'))
    Run_filelist = sorted(Run_filelist)
    # Pick the first file from the list (assuming only one file per directory)
    RunFile = Run_filelist[0]
   
    # Read in the dataset from the file
    Run_ds = xr.open_dataset(RunFile)
    print('Read in files %s' % (RunFile))
    
    # Select the variable data for the specified vertical level from the dataset
    Run_DayLevi_Vari = Run_ds.isel(lev=lev_idx)
    
    if varname!='None':
        Run_DayLevi_Vari = Run_DayLevi_Vari[varname]
    
    print('Run - Vertical Level: %s' % (str(int(Run_DayLevi_Vari.lev.values)) + ' hPa'))
    
    return Run_DayLevi_Vari


########### 
import xarray as xr
import numpy as np
import glob

def comp_2RunResults_h1_diff(Run1Path, Run2Path, datei, lev_idx, varname):
    """
    Compare two datasets for a given variable at a specific time and vertical level.
    Run1 - Run2

    Parameters:
    Run1Path (str): Path to the directory containing the first set of files.
    Run2Path (str): Path to the directory containing the second set of files.
    histfreq (str): Keyword used to identify the files (e.g., 'hist', 'day', 'hourly', etc.).
    time_idx (int): Index of the desired time slice.
    lev_idx (int): Index of the desired vertical level (e.g., pressure level).
    varname (str): Name of the variable to compare.

    Returns:
    dict: A dictionary containing the following:
          - 'ABSdiff': Absolute difference between the two datasets.
          - 'RELAdiff': Relative difference as a percentage relative to the first dataset.
          - 'unit': Units of the variable.
          - 'longname': Long name/description of the variable.
          - 'Timei': Time value of the selected time index in the first dataset.
    """
    # for daily average file
    histfreq = 'h1'
    datei_str = datei+'T00:00:00.000000000'
    
    # Get a list of file paths in the specified directories with the given keyword
    Run1_filelist = glob.glob(os.path.join(Run1Path, f'*{histfreq}*'))
    Run2_filelist = glob.glob(os.path.join(Run2Path, f'*{histfreq}*'))

    # Pick the first file from each list (assuming only one file per directory)
    Run1File = Run1_filelist[0]
    Run2File = Run2_filelist[0]
    
    # Read in the datasets from the files
    Run1_ds = xr.open_dataset(Run1File)
    Run2_ds = xr.open_dataset(Run2File)
    print('Read in files %s and %s' % (Run1File, Run2File))
    
    # Select the variable data for the specified time and vertical level from each dataset
    Run1_DayLevi_Vari = Run1_ds.isel(lev=lev_idx).sel(time=datei_str)[varname]
    print('R1 - Time: %s; Vertical Level: %s' % (str(Run1_DayLevi_Vari.time.values),
                                                 str(int(Run1_DayLevi_Vari.lev.values)) + ' hPa'))
    
    Run2_DayLevi_Vari = Run2_ds.isel(lev=lev_idx).sel(time=datei_str)[varname]
    print('R2 - Time: %s; Vertical Level: %s' % (str(Run2_DayLevi_Vari.time.values),
                                                 str(int(Run2_DayLevi_Vari.lev.values)) + ' hPa'))
    
    # Calculate the difference between the two datasets
    difference_ds = Run2_DayLevi_Vari-Run1_DayLevi_Vari
    
    # Absolute difference
    ABSdiff = difference_ds[:].values
    
    # Relative difference to Run1
    RELAdiff = 100 * np.divide(difference_ds[:].values, Run1_DayLevi_Vari[:].values)
    
    # Get the unit and long name of the variable
    unit = Run1_DayLevi_Vari.units
    longname = Run1_DayLevi_Vari.long_name
    
    # Return the results in a dictionary
    return {'ABSdiff': ABSdiff, 'RELAdiff': RELAdiff, 'unit': unit,
            'longname': longname, 'Timei': str(Run1_DayLevi_Vari.time.values)}


########### 
import xarray as xr
import numpy as np
import glob

def comp_2RunResults_diff(Run1Path, Run2Path, histfreq, time_idx, lev_idx, varname):
    """
    Compare two datasets for a given variable at a specific time and vertical level.

    Parameters:
    Run1Path (str): Path to the directory containing the first set of files.
    Run2Path (str): Path to the directory containing the second set of files.
    histfreq (str): Keyword used to identify the files (e.g., 'hist', 'day', 'hourly', etc.).
    time_idx (int): Index of the desired time slice.
    lev_idx (int): Index of the desired vertical level (e.g., pressure level).
    varname (str): Name of the variable to compare.

    Returns:
    dict: A dictionary containing the following:
          - 'ABSdiff': Absolute difference between the two datasets.
          - 'RELAdiff': Relative difference as a percentage relative to the first dataset.
          - 'unit': Units of the variable.
          - 'longname': Long name/description of the variable.
          - 'Timei': Time value of the selected time index in the first dataset.
    """
    
    # Get a list of file paths in the specified directories with the given keyword
    Run1_filelist = glob.glob(os.path.join(Run1Path, f'*{histfreq}*'))
    Run2_filelist = glob.glob(os.path.join(Run2Path, f'*{histfreq}*'))

    # Pick the first file from each list (assuming only one file per directory)
    Run1File = Run1_filelist[0]
    Run2File = Run2_filelist[0]
    
    # Read in the datasets from the files
    Run1_ds = xr.open_dataset(Run1File)
    Run2_ds = xr.open_dataset(Run2File)
    print('Read in files %s and %s' % (Run1File, Run2File))
    
    # Select the variable data for the specified time and vertical level from each dataset
    Run1_DayLevi_Vari = Run1_ds.isel(time=time_idx, lev=lev_idx)[varname]
    print('R1 - Time: %s; Vertical Level: %s' % (str(Run1_DayLevi_Vari.time.values),
                                                 str(int(Run1_DayLevi_Vari.lev.values)) + ' hPa'))
    
    Run2_DayLevi_Vari = Run2_ds.isel(time=time_idx, lev=lev_idx)[varname]
    print('R2 - Time: %s; Vertical Level: %s' % (str(Run2_DayLevi_Vari.time.values),
                                                 str(int(Run2_DayLevi_Vari.lev.values)) + ' hPa'))
    
    # Calculate the difference between the two datasets
    difference_ds = Run2_DayLevi_Vari-Run1_DayLevi_Vari
    
    # Absolute difference
    ABSdiff = difference_ds[:].values
    
    # Relative difference to Run1
    RELAdiff = 100 * np.divide(difference_ds[:].values, Run1_DayLevi_Vari[:].values)
    
    # Get the unit and long name of the variable
    unit = Run1_DayLevi_Vari.units
    longname = Run1_DayLevi_Vari.long_name
    
    # Return the results in a dictionary
    return {'ABSdiff': ABSdiff, 'RELAdiff': RELAdiff, 'unit': unit,
            'longname': longname, 'Timei': str(Run1_DayLevi_Vari.time.values)}


########### 


########### 



########### 


########### 



########### 