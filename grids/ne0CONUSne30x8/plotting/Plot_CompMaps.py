# This script plots the comparison maps for Svante and Cheyenne runs

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../functions/'))

import numpy as np # for array manipulation and basic scientific calculation
import xarray as xr # To read NetCDF files
import matplotlib.pyplot as plt # Core library for plotting
import matplotlib.cm as cm # To use different colormaps
from Plot_2D import Plot_2D # To draw a map
import cartopy.crs as ccrs # For map projection
import seaborn as sns # boxplot
import pandas as pd

# save PDF
from matplotlib.backends.backend_pdf import PdfPages

# use time new roman for all text
plt.rcParams["font.family"] = "Times New Roman"

import warnings

# Ignore a specific warning
warnings.filterwarnings('ignore', message='Some warning message')
# Ignore the FutureWarning caused by iteritems()
warnings.filterwarnings("ignore", category=FutureWarning, message=".*iteritems.*")

# ============================================================
# CONFIGURATION - every path comes from config/paths.py, the single
# source of truth. Override cluster locations with the MUSICA_ENV_*
# environment variables documented there. Do not hard-code paths here.
# ============================================================
import sys as _sys, pathlib as _pathlib
_ROOT = next(_p for _p in _pathlib.Path(__file__).resolve().parents
             if (_p / 'config' / 'paths.py').exists())
_sys.path.insert(0, str(_ROOT))
import config  # noqa: F401  - also puts functions/ on sys.path
from config.paths import (
    ARCHIVE,
    CHEYENNE_ARCHIVE,
    CESM_FIGURE_DIR,
    ensure_dir,
)

cheyenne_archive = str(CHEYENNE_ARCHIVE) + '/'
svante_archive = str(ARCHIVE) + '/'
figure_diri = str(ensure_dir(CESM_FIGURE_DIR / 'CompMaps')) + '/'

casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.base.001'

# ============================================================

# Read SCRIP file that has grid information needed to plot values on a map
SCRIP_CONUS = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../grid_files/ne0CONUS_ne30x8_np4_SCRIP.nc')
SCRIP_ne30 = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../grid_files/ne30np4_091226_pentagons.nc')

cheyenne_filepath = cheyenne_archive+casename+'/atm/hist/'+casename+'.cam.h1.2013-08-01-00000.nc'
cheyenne_ds_CONUS = xr.open_dataset( cheyenne_filepath ) 

svante_filepath = svante_archive+casename+'/atm/hist/'+casename+'.cam.h1.2013-08-01-00000.nc'
svante_ds_CONUS = xr.open_dataset( svante_filepath ) 


###################
"""Absolute difference - A series of difference plots for all days into one PDF"""
varlist = ['PS','Q','T','CLOUD','PM25_SRF']
Region = 'CONUS'
DiffType = 'Absolute Difference' # Absolute(ABS) or Relative (RELA)

PDFOut_name = 'CompSvanteChey_ABSDiff_Map'+Region+'.pdf'

# need to generate n figures
nfig = len(cheyenne_ds_CONUS.time.values)*len(varlist)

# create a list of plots to save
plots = []
for figi in range(nfig):
    """Plot"""
    # for each variable
    for varname in varlist:
        # loop through all time
        for timei in cheyenne_ds_CONUS.time.values:
            print(varname,timei)
            Datei = str(timei)[:10]
            cheyenne_selda = cheyenne_ds_CONUS.sel(time=timei)
            svante_selda = svante_ds_CONUS.sel(time=timei)

            ### Get difference array
            if cheyenne_selda[varname].ndim==2:
                # -1 for surface, : for global (all columns) 
                # absolute difference
                ABSdiff_Dayi = svante_selda[varname][-1,:].values-cheyenne_selda[varname][-1,:].values
                # relative difference to cheyenne
                RELAdiff_Dayi = 100*np.divide((svante_selda[varname][-1,:].values-cheyenne_selda[varname][-1,:].values),cheyenne_selda[varname][-1,:].values)
            elif cheyenne_selda[varname].ndim==1:
                # already at the surface 
                # absolute difference
                ABSdiff_Dayi = svante_selda[varname][:].values-cheyenne_selda[varname][:].values
                # relative difference to cheyenne
                RELAdiff_Dayi = 100*np.divide((svante_selda[varname][:].values-cheyenne_selda[varname][:].values),cheyenne_selda[varname][:].values)

            # set up for absolute versus relative difference
            if DiffType=='Absolute Difference':
                Plot_ar = ABSdiff_Dayi
                Plot_unit = cheyenne_selda[varname].units
            elif DiffType=='Relative Difference':
                Plot_ar = RELAdiff_Dayi
                Plot_unit = '%'
                # rangeMax = 100

            # need to set the rangeMax, but avoid giving all zeros if the difference is small
            if (np.nanmax(Plot_ar))<1:
                rangeMax = 1
            else:
                rangeMax = np.nanmax(Plot_ar)

            ### Which map
            fig = plt.figure( figsize=(16,14) ) 
            # - ne30x8 regional refinement over CONUS|
            ax1 = fig.add_subplot(1,1,1,projection=ccrs.PlateCarree())
            if Region=="CONUS":
                Plot_2D( Plot_ar, scrip_file=SCRIP_CONUS, ax=ax1,
                        cmin=-rangeMax, cmax=rangeMax, 
                        unit=Plot_unit,
                        cmap=cm.bwr, state=True, lon_range=[-140,-50], lat_range=[15,60],
                          grid_line=False, grid_line_lw=0.15 ) # cmap=cm.hot_r, 
            elif Region=="Global":
                Plot_2D( Plot_ar, scrip_file=SCRIP_CONUS, ax=ax1,
                        cmin=-rangeMax, cmax=rangeMax, 
                        unit=Plot_unit,
                        cmap=cm.bwr, state=True,
                          grid_line=False, grid_line_lw=0.15 ) # cmap=cm.hot_r, 
            # finish plotting
            plt.title(DiffType+' in '+cheyenne_selda[varname].long_name+' \n'+Datei+' (cheyenne - svante)', fontsize=16, y=1.02);
            plots.append(fig)
            
# create a new PDF file and save the plots to separate pages
with PdfPages(PDFOut_name) as pdf:
    for plot in plots:
        pdf.savefig(plot)

# close all figures
plt.close('all')
print('save to: '+PDFOut_name)

###################
"""Relative difference - A series of difference plots for all days into one PDF"""
varlist = ['PS','Q','T','CLOUD','PM25_SRF']
Region = 'CONUS'
DiffType = 'Relative Difference' # Absolute(ABS) or Relative (RELA)

PDFOut_name = 'CompSvanteChey_RELADiff_Map'+Region+'.pdf'

# need to generate n figures
nfig = len(cheyenne_ds_CONUS.time.values)*len(varlist)

# create a list of plots to save
plots = []
for figi in range(nfig):
    """Plot"""
    # for each variable
    for varname in varlist:
        # loop through all time
        for timei in cheyenne_ds_CONUS.time.values:
            print(varname,timei)
            Datei = str(timei)[:10]
            cheyenne_selda = cheyenne_ds_CONUS.sel(time=timei)
            svante_selda = svante_ds_CONUS.sel(time=timei)

            ### Get difference array
            if cheyenne_selda[varname].ndim==2:
                # -1 for surface, : for global (all columns) 
                # absolute difference
                ABSdiff_Dayi = svante_selda[varname][-1,:].values-cheyenne_selda[varname][-1,:].values
                # relative difference to cheyenne
                RELAdiff_Dayi = 100*np.divide((svante_selda[varname][-1,:].values-cheyenne_selda[varname][-1,:].values),cheyenne_selda[varname][-1,:].values)
            elif cheyenne_selda[varname].ndim==1:
                # already at the surface 
                # absolute difference
                ABSdiff_Dayi = svante_selda[varname][:].values-cheyenne_selda[varname][:].values
                # relative difference to cheyenne
                RELAdiff_Dayi = 100*np.divide((svante_selda[varname][:].values-cheyenne_selda[varname][:].values),cheyenne_selda[varname][:].values)

            # set up for absolute versus relative difference
            if DiffType=='Absolute Difference':
                Plot_ar = ABSdiff_Dayi
                Plot_unit = cheyenne_selda[varname].units
            elif DiffType=='Relative Difference':
                Plot_ar = RELAdiff_Dayi
                Plot_unit = '%'
                # rangeMax = 100

            # need to set the rangeMax, but avoid giving all zeros if the difference is small
            if (np.nanmax(Plot_ar))<1:
                rangeMax = 1
            else:
                rangeMax = np.nanmax(Plot_ar)

            ### Which map
            fig = plt.figure( figsize=(16,14) ) 
            # - ne30x8 regional refinement over CONUS|
            ax1 = fig.add_subplot(1,1,1,projection=ccrs.PlateCarree())
            if Region=="CONUS":
                Plot_2D( Plot_ar, scrip_file=SCRIP_CONUS, ax=ax1,
                        cmin=-rangeMax, cmax=rangeMax, 
                        unit=Plot_unit,
                        cmap=cm.bwr, state=True, lon_range=[-140,-50], lat_range=[15,60],
                          grid_line=False, grid_line_lw=0.15 ) # cmap=cm.hot_r, 
            elif Region=="Global":
                Plot_2D( Plot_ar, scrip_file=SCRIP_CONUS, ax=ax1,
                        cmin=-rangeMax, cmax=rangeMax, 
                        unit=Plot_unit,
                        cmap=cm.bwr, state=True,
                          grid_line=False, grid_line_lw=0.15 ) # cmap=cm.hot_r, 
            # finish plotting
            plt.title(DiffType+' in '+cheyenne_selda[varname].long_name+' \n'+Datei+' (cheyenne - svante)', fontsize=16, y=1.02);
            plots.append(fig)
            
# create a new PDF file and save the plots to separate pages
with PdfPages(PDFOut_name) as pdf:
    for plot in plots:
        pdf.savefig(plot)

# close all figures
plt.close('all')
print('save to: '+PDFOut_name)