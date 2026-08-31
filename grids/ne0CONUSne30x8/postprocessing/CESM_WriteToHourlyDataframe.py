#===============
# This script writes all CESM variables/Given varaibles at one site for given dates in EDT to a CSV file

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
    AQS_DIR,
    WRFCMAQ_LISTOS_DIR,
    CMAQ_ACONC_DIR,
    SITE_HOURLY_JJA2018_DIR,
    ensure_dir,
)

AQS_diri = str(AQS_DIR) + '/'
CESM_diri = str(WRFCMAQ_LISTOS_DIR) + '/'
SiteCSV_diri = str(ensure_dir(SITE_HOURLY_JJA2018_DIR)) + '/'

# ============================================================

# Need to replace!
Datelist = ['20180601','20180602','20180603','20180604','20180605','20180625','20180630',
           '20180701','20180702','20180703','20180719','20180720']

#===============
# Functions from the function dir
import sys, os, glob
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../functions/'))

# Legacy imports from the original functions directory; these may not be available in all environments
try:
    from external_packages import *
except ImportError:
    pass  # optional legacy import — install or provide external_packages.py if needed
try:
    from plot_settings import *
except ImportError:
    pass  # optional legacy import — install or provide plot_settings.py if needed
try:
    from statistics_func import *
except ImportError:
    pass  # optional legacy import — install or provide statistics_func.py if needed

# additional library for writing MATLAB files
import scipy.io as sio

#===============
def get_current_next_days(date_string):
    from datetime import datetime, timedelta
    """Get the nearest three days centered on a given date"""
    # Convert input date to datetime object
    date_obj = datetime.strptime(date_string, '%Y%m%d')

    # Calculate the dates for the previous day, current day, and next day
    previous_day = (date_obj - timedelta(days=1)).strftime('%Y%m%d')
    current_day = date_obj.strftime('%Y%m%d')
    next_day = (date_obj + timedelta(days=1)).strftime('%Y%m%d')

    # Combine the dates into two different lists with different formats
    given_days_format_1 = [current_day, next_day]
    given_days_format_2 = [
        datetime.strptime(day, '%Y%m%d').strftime('%Y-%m-%d')
        for day in given_days_format_1
    ]

    # Return the resulting lists
    return given_days_format_1, given_days_format_2

# handle datetime
def dt64_to_datetime(dft6):
    """Converts dt64 to datetime"""
    from datetime import datetime
    return datetime.utcfromtimestamp(dft6.tolist()/1e9)


#===============
#===============
#===============

# directory locations (CESM_diri configured in USER CONFIGURATION block above)
ACONC_diri = str(CMAQ_ACONC_DIR) + '/'

lat2d = np.load(CESM_diri+'CMAQ_lat2d.npy')
lon2d = np.load(CESM_diri+'CMAQ_lon2d.npy')

# match the keys
sitedf = pd.read_csv(AQS_diri+"SiteSetting_Pandoradf_2018_CT_LISTOS_summer.csv")
# sitedf
Pandora_colocate_Monitor_dic = dict(zip(sitedf.SitesLocs.values, sitedf.SiteAbbr.values))
Pandora_SiteLoc_dic = dict(zip(sitedf.SiteAbbr.values, sitedf.SitesLocs.values))
AQSID_Pandora_SitenameDic = dict(zip(sitedf.SitesLocs.values, sitedf.Pandora_SiteNames.values))
LISTOS_Pandora_SitenameDic = dict(zip(sitedf.SiteAbbr.values, sitedf.Pandora_SiteNames.values))

# lat & lon
Pandora_SiteLat_dic = dict(zip(sitedf.SiteAbbr.values, sitedf.Latitude.values))
Pandora_SiteLon_dic = dict(zip(sitedf.SiteAbbr.values, sitedf.Longitude.values))

## get ROW and COL index
sitedf = pd.read_csv(CESM_diri+"CMAQIdx_SiteSetting_Pandoradf_2018_CT_LISTOS_summer.csv")

ROWidx_dic = dict(zip(sitedf.SiteAbbr.values, sitedf.CMAQ_ROWidx.values))
COLidx_dic = dict(zip(sitedf.SiteAbbr.values, sitedf.CMAQ_COLidx.values))

# read in one example file


def aconc_path(datei):
    """Resolve the CMAQ ACONC file for one YYYYMMDD date.

    The delivered file names carry a trailing per-recipient suffix, so match on
    the date by glob rather than constructing an exact file name.
    """
    hits = sorted(glob.glob(os.path.join(ACONC_diri, f'CCTM_ACONC_*_{datei}*.nc')))
    if not hits:
        raise FileNotFoundError(
            f'No CMAQ ACONC file for {datei} in {ACONC_diri}')
    return hits[0]


# read in data
testCMAQ_DA = xr.open_dataset(aconc_path('20180601'))#.sel(TSTEP=slice(0,24))
TSTEPVals = testCMAQ_DA.TSTEP.values
# print(TSTEPVals)

variable_names = list(testCMAQ_DA.variables.keys())
# remove the string "TFLAG" from the list
variable_names = [var_name for var_name in testCMAQ_DA.variables if var_name != 'TFLAG']

def ACONC_getvari_siteDA_UTC(ACONC_diri,datei,Sitei,varname):
    """Get data for selected day
        Prefix:
            "ACONC"
    """    
    da = xr.open_dataset(aconc_path(datei)).sel(TSTEP=slice(0,24))
    
    ROWidx = ROWidx_dic[Sitei]
    COLidx = COLidx_dic[Sitei]

    # Site-selected pixel dataset
    site_da = da.sel(COL=COLidx,ROW=ROWidx,LAY=0)  
    
    # for the given variable
    site_var_da = site_da[varname]    
    site_var = site_var_da.values.flatten()
    var_unit = site_var_da.units
    var_unit = var_unit.replace(" ", "") # get rid of spaces
    
    return {'site_var':site_var,'var_unit':var_unit}

###############
Prefix = "CMAQSurfACONCallChemicals"
# SiteCSV_diri configured in USER CONFIGURATION block above

# all variables in one dataframe Per site
for Sitei in sitedf.SiteAbbr.values:
    print(Sitei)
    CVSfileout_name = Sitei+"_"+Prefix+"_SelectedDays_24hrEDT.csv"
    Sitei_alldays_allvars_df = pd.DataFrame()
    
    for dateidx in range(len(Datelist)):
        Filefmt_Datei = Datelist[dateidx]
        date_obj = datetime.strptime(Filefmt_Datei, '%Y%m%d')
        Datei = date_obj.strftime('%Y-%m-%d')
        
        # for this day
        Sitei_datei_allvars_df = pd.DataFrame()

        for varidx in range(len(variable_names)):
            
            varname = variable_names[varidx]
            # get the data for each day-24hr
            site_var_dic = ACONC_getvari_siteDA_UTC(ACONC_diri,Filefmt_Datei,Sitei,varname)
            site_var = site_var_dic['site_var']
            var_unit = site_var_dic['var_unit']

            # Convert to EDT and then get the date time as separate str for the selected days
            datetimes = []
            DateLocals = []
            TimeLocals = []

            from datetime import datetime
            from pytz import timezone

            for TSTEPi in TSTEPVals:
                date_str = str(Datei)+" "+str(TSTEPi).zfill(2)+":00:00 UTC+0000"
                in_utc = datetime.strptime(date_str, "%Y-%m-%d %H:%M:%S %Z%z") #dt64_to_datetime(dt64i)
                in_edt = in_utc.astimezone(timezone('US/Eastern'))
                datetimes.append(in_edt.strftime("%Y-%m-%d %H:%M:%S"))
                DateLocals.append(in_edt.strftime("%Y-%m-%d"))
                TimeLocals.append(in_edt.strftime("%H:%M"))    

            # write to dataframe
            VarColName = Prefix+"_"+varname+"_"+var_unit
            Sitei_varidf_24hr_df = pd.DataFrame(list(zip([Sitei]*len(datetimes),datetimes,DateLocals,TimeLocals,site_var)), 
                                              columns=['site','datetime','Date Local','Time Local',VarColName])
            
            # merge for different variables (need to begin logging with the merged columns)
            if varidx==0:
                Sitei_datei_allvars_df = Sitei_varidf_24hr_df
            else:
                Sitei_datei_allvars_df = pd.merge(Sitei_datei_allvars_df, Sitei_varidf_24hr_df, how='left', 
                                         on=['site','datetime',"Date Local","Time Local"])
            
        # append for each day to the end
        Sitei_alldays_allvars_df = Sitei_alldays_allvars_df.append(Sitei_datei_allvars_df, ignore_index = True)

    # save to dataframe
    Sitei_alldays_allvars_df.to_csv(SiteCSV_diri+CVSfileout_name, index=False)
    print("Save to:",CVSfileout_name)
