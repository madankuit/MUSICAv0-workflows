"""
config/paths.py — the single source of truth for every path in this repository.

Rules
-----
1. No script anywhere in this repo may hard-code a path. Scripts import from here.
2. Files that ship WITH the repo (SCRIP grids, masks, target-grid descriptions,
   shared functions) are addressed **repo-relatively** and therefore work on any
   machine, in any checkout location, with no configuration at all.
3. Files that live on a cluster (model archive, satellite data, emissions,
   processed output) are built from the *runtime* user name, so this file
   contains no personal name and is correct for whoever runs it.

Where new cluster output belongs
--------------------------------
Cluster paths fall into two groups, and the distinction matters when adding one:

* **This project's own data** -> under ``MUSICA_PROJECT_ROOT``
  (``$DATA_ROOT/MUSICAv0-workflows/{data,figures}/``), mirroring the layout the
  TEMPO_OTR_2024O3Events project uses. **Put new MUSICA-specific outputs here.**
* **Shared, cross-project data** -> the data-type folders (``CESM22/``,
  ``Datasets/``, ``ncar_copies/``): the model archive, raw and regridded
  satellite products, emission inventories. These are inputs other projects also
  read, so they are not project-scoped.

Some existing outputs predate this convention and deliberately stay where they
are, because moving them would break paths for no scientific gain --
``BGO3_ROOT`` (under ``ProcessedData/``),
``TROP_VCD_ROOT``/``TOTAL_VCD_ROOT`` and ``REGRIDDED_2018_ROOT`` (under
``CESM22/``). Follow the rule above for anything new rather than matching them.

Overriding the cluster locations
--------------------------------
Everything under "MACHINE-SPECIFIC ROOTS" is derived from two roots that can be
overridden with environment variables — set these if your data does not live
under your own user name::

    export MUSICA_ENV_DATA_ROOT=/net/fs09/d0/<somebody>     # large data volume
    export MUSICA_ENV_HOME_ROOT=/home/<somebody>            # home / scripts volume

Individual leaf paths can be overridden the same way; see `_env_path` uses below.

Usage
-----
Any script in this repo reaches the config with a depth-independent bootstrap
(no absolute paths, works from any subdirectory)::

    import sys, pathlib
    _root = next(p for p in pathlib.Path(__file__).resolve().parents
                 if (p / "config" / "paths.py").exists())
    sys.path.insert(0, str(_root))
    import config                      # also puts functions/ on sys.path
    from config.paths import ARCHIVE, SCRIP_NE30NP4

Catalogue entries
-----------------
A handful of constants are declared here but not read by any committed script or
notebook: ``CAMS_V51_ORIG_DIR``, ``CAMS_V62_ORIG_DIR``,
``CAMS_V62_NE30NP4_NONAN_DIR``, ``CAMS_NEI2017_MERGED_ROOT``, ``NEI2017_01_DIR``,
``NEI2022V2_01_DIR``, ``TROPOMI_L2_HCHO_DIR``, ``MUSICA_PROJECT_FIGURES``.
They are **not dead code** — they are inherited from the reference file this
module replaced, they record where those inputs actually live, and they all
resolve on the reference cluster (verified via ``check_paths.py``). Keep them
unless the underlying data genuinely moves or disappears.

MODIFICATION HISTORY:
    VERSION 1.0
    - Initial version; replaces svante_MUSICA_paths.py. All paths consolidated
      here, repo-shipped data addressed repo-relatively, user name removed.
"""

import os
import getpass
from pathlib import Path


def _env_path(var, default):
    """Return Path(os.environ[var]) if set, else `default`."""
    val = os.environ.get(var)
    return Path(val) if val else Path(default)


# ============================================================
# 1. REPO-RELATIVE PATHS  —  data that ships with this repository
#    These need no configuration and are valid on every machine.
# ============================================================

REPO_ROOT = Path(__file__).resolve().parents[1]

FUNCTIONS_DIR = REPO_ROOT / 'functions'
GRIDS_DIR = REPO_ROOT / 'grids'

# --- ne30np4: standard global SE grid (~111 km, 48,602 ncol) ---
NE30NP4_GRID_DIR = GRIDS_DIR / 'ne30np4' / 'grid_files'
SCRIP_NE30NP4 = NE30NP4_GRID_DIR / 'ne30np4_091226_pentagons.nc'
MASK_NE30NP4_CONUS_80KM = NE30NP4_GRID_DIR / 'ne30np4_091226_pentagons_CONUSlandMaskedFalse_80kmBuffer.nc'
MASK_NE30NP4_LOWER48_50KM = NE30NP4_GRID_DIR / 'ne30np4_091226_pentagons_Lower48StatesCoastal50kmMaskedFalse.nc'

# --- ne0CONUSne30x8: variable-resolution SE grid (~14 km CONUS, 174,098 ncol) ---
NE0CONUS_GRID_DIR = GRIDS_DIR / 'ne0CONUSne30x8' / 'grid_files'
SCRIP_NE0CONUSNE30X8 = NE0CONUS_GRID_DIR / 'ne0CONUS_ne30x8_np4_SCRIP.nc'
MASK_NE0CONUS_CONUS_80KM = NE0CONUS_GRID_DIR / 'ne0CONUSne30x8_np4_CONUSlandMaskedFalse_80kmBuffer.nc'
MASK_NE0CONUS_US_REGIONS = NE0CONUS_GRID_DIR / 'US_region_masks_ne0CONUSne30x8_ncol174098.nc'
MASK_NE0CONUS_CANADA_PROVINCES = NE0CONUS_GRID_DIR / 'mask_CanadaProvinces_ne0CONUS_ne30x8.nc'
MASK_NE0CONUS_CANADIAN_FIRE = NE0CONUS_GRID_DIR / 'ne0CONUSne30x8_Canadianfire_region_masks.nc'
MASK_F09_CONUS_80KM = NE0CONUS_GRID_DIR / 'f09_CONUSlandMaskedFalse_80kmBuffer.nc'

# --- regular lat/lon destination grids used as ESMF regrid targets ---
TARGET_GRID_DIR = GRIDS_DIR / 'target_grids' / 'grid_files'
FV_GRIDINFO_1X1 = TARGET_GRID_DIR / 'FV1x1grid_info_c20241105.nc'
FV_GRIDINFO_015 = TARGET_GRID_DIR / 'FV_gridinfo_0.15_c20231204.nc'
FV_GRIDINFO_01_CAMS = TARGET_GRID_DIR / 'FV_gridinfo_CAMS_c20210219.nc'


# ============================================================
# 2. MACHINE-SPECIFIC ROOTS  —  cluster data, not in this repository
#    Derived from the runtime user; override with the env vars above.
# ============================================================

USER = os.environ.get('USER') or os.environ.get('LOGNAME') or getpass.getuser()

DATA_ROOT = _env_path('MUSICA_ENV_DATA_ROOT', f'/net/fs09/d0/{USER}')
HOME_ROOT = _env_path('MUSICA_ENV_HOME_ROOT', Path.home())

# Project folder for this repository's own data and figures, mirroring the
# TEMPO_OTR_2024O3Events layout (data/ + figures/). New MUSICA-specific outputs
# belong here rather than in the shared data-type folders below.
MUSICA_PROJECT_ROOT = _env_path('MUSICA_ENV_PROJECT_ROOT',
                                DATA_ROOT / 'MUSICAv0-workflows')
MUSICA_PROJECT_DATA = MUSICA_PROJECT_ROOT / 'data'
MUSICA_PROJECT_FIGURES = MUSICA_PROJECT_ROOT / 'figures'

CESM22_ROOT = DATA_ROOT / 'CESM22'
DATASETS_ROOT = DATA_ROOT / 'Datasets'
PROCESSED_DATA_ROOT = DATA_ROOT / 'ProcessedData'
NCAR_COPIES_ROOT = DATA_ROOT / 'ncar_copies'
FIGURES_ROOT = DATA_ROOT / 'Figures'
CESM_FIGURE_DIR = FIGURES_ROOT / 'CESM_analysis'


# ============================================================
# 3. ESMF REGRIDDING WEIGHTS  —  large derived files, kept on the cluster
#    Regenerate with the gen_*_weights.py scripts if missing.
# ============================================================

GRIDS_EXTERNAL_DIR = _env_path('MUSICA_ENV_GRIDS_DIR', CESM22_ROOT / 'grids')

# ne0CONUSne30x8 -> regular lat/lon
WEIGHTS_NE0CONUS_TO_015 = GRIDS_EXTERNAL_DIR / 'ne0CONUSne30x8_ESMFmap_0.15x0.15_cubit_conserve_cams_c20231204.nc'
WEIGHTS_NE0CONUS_TO_01 = GRIDS_EXTERNAL_DIR / 'ESMFmap_0.1x0.1_ne0CONUSne30x8_cubit_conserve_cams.nc'
# ne30np4 <-> regular lat/lon
WEIGHTS_FV09X125_TO_NE30 = GRIDS_EXTERNAL_DIR / 'ESMFmap_0.9x1.25_ne30np4_cubit_conserve_cams.nc'
WEIGHTS_NE30_TO_1X1 = _env_path(
    'MUSICA_ENV_WEIGHTS_NE30_TO_1X1',
    GRIDS_EXTERNAL_DIR / 'ESMFmap_ne30np4_TO_1x1_conserve_c20260708.nc')

# Alternative ne30np4 SCRIP (5 MB, not shipped in-repo; optional)
SCRIP_NE30NP4_ALT = GRIDS_EXTERNAL_DIR / 'ne30np4_grid_c20241105.nc'


# ============================================================
# 4. CESM / MUSICA MODEL OUTPUT
# ============================================================

ARCHIVE = CESM22_ROOT / 'archive'
# Archive of runs carried over from NCAR (Cheyenne/Derecho), used for run-vs-run maps
CHEYENNE_ARCHIVE = _env_path('MUSICA_ENV_CHEYENNE_ARCHIVE', NCAR_COPIES_ROOT / 'archive')

PROCESSED_OUTPUT_DIR = CESM22_ROOT / 'processed_output'
JULY2018_SURFH2_MERGED_DIR = PROCESSED_OUTPUT_DIR / 'July2018_surfh2_merged'

# Vertical column densities on the native ncol grid; append the casename
TROP_VCD_ROOT = CESM22_ROOT / 'Calculated_MUSICA_VCD' / 'TroposphericVCD'
TOTAL_VCD_ROOT = CESM22_ROOT / 'Calculated_MUSICA_VCD' / 'TotalVCD'

# External meteorology file for cases whose h2 stream did not save PS
MET_FILE_JULY2018 = (TROP_VCD_ROOT /
                     'met_f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1base01.nc')

# Regridded to 0.15 deg lat/lon (July 2018, TROPOMI overpass time)
REGRIDDED_2018_ROOT = (CESM22_ROOT / 'Regridded_MUSICA_Output' /
                       '2018_1330LT_TROPOMIcomp' / 'MassConserve_latlon015_MUSICAoutput')

# Per-case subdirectory names under REGRIDDED_2018_ROOT / <casename> /
REGRID_SUBDIR_H2 = 'h2'
REGRID_SUBDIR_BYVAR = 'h2_ByVar'
REGRID_SUBDIR_VERT_BYVAR_DATE = 'h2_VertiVarsByVarDate'
REGRID_SUBDIR_TROPVCD = 'TropVCD_approx1330LT'
REGRID_SUBDIR_TOTALVCD = 'TotalVCD_approx1330LT'


def case_hist_dir(casename, archive=None):
    """Directory holding the raw CAM history files for `casename`."""
    return Path(archive or ARCHIVE) / casename / 'atm' / 'hist'


def case_regrid_dir(casename, subdir=REGRID_SUBDIR_H2, root=None):
    """Directory holding regridded output for `casename` (see REGRID_SUBDIR_*)."""
    return Path(root or REGRIDDED_2018_ROOT) / casename / subdir


def case_trop_vcd_dir(casename, root=None):
    """Directory holding native-grid tropospheric VCDs for `casename`."""
    return Path(root or TROP_VCD_ROOT) / casename


def case_total_vcd_dir(casename, root=None):
    """Directory holding native-grid total VCDs for `casename`."""
    return Path(root or TOTAL_VCD_ROOT) / casename


# ============================================================
# 5. TROPOMI SATELLITE DATA
# ============================================================

TROPOMI_L2_NO2_DIR = DATASETS_ROOT / 'S5P_L2__NO2____HiR_2' / '2018'
TROPOMI_L2_HCHO_DIR = DATASETS_ROOT / 'S5P_L2__HCHO___HiR_2' / '2018'
TROPOMI_L2_CO_DIR = DATASETS_ROOT / 'S5P_L2__CO_____HiR_2' / '2018'

# Regridded to 0.15 deg, matching the conservative MUSICA regrid
TROPOMI_015_NO2_DIR = DATASETS_ROOT / 'S5P_L2__NO2____HiR_2' / 'regrid_2018_MUSICA015'
TROPOMI_015_HCHO_DIR = DATASETS_ROOT / 'S5P_L2__HCHO___HiR_2' / 'regrid_2018_MUSICA015'
TROPOMI_015_CO_DIR = DATASETS_ROOT / 'S5P_L2__CO_____HiR_2' / 'regrid_2018_MUSICA015'

# Averaging kernels and AMFs regridded to 0.1 deg global
_TROPOMI_L3_ROOT = DATASETS_ROOT / 'TROPOMI_RegridL2_ToL3' / 'Global'
TROPOMI_01_NO2_AK_DIR = _TROPOMI_L3_ROOT / 'TROPOMI_NO2AK_01deg'
TROPOMI_01_HCHO_AK_DIR = _TROPOMI_L3_ROOT / 'TROPOMI_HCHOAK_01deg'
TROPOMI_01_NO2_DIR = _TROPOMI_L3_ROOT / 'TROPOMI_NO2_01deg'
# Required to convert the stored total-column NO2 AK to a tropospheric AK
# (AK_trop = AK_total * AMF_total / AMF_trop); see functions/README.md.
#
# These two do NOT exist yet — the only outstanding gap reported by check_paths.py.
# NOT a data-acquisition problem: air_mass_factor_total and
# air_mass_factor_troposphere are ordinary TROPOMI L2 variables sitting in the
# same /PRODUCT/ group as averaging_kernel. The L2->L3 regrid that produced the AK
# files simply carried the AK only, so closing this means re-running that regrid
# with the two AMF variables added. (HCHO and CO are unaffected.) They are deliberately kept HERE, beside TROPOMI_NO2AK_01deg and
# the other regridded TROPOMI L3 products, rather than under
# MUSICA_PROJECT_ROOT: the vertical-regrid code reads AK and AMF together per
# pixel, so separating them makes the pairing easy to break, and these are
# shared satellite inputs rather than MUSICA outputs.
TROPOMI_01_NO2_AMF_TOTAL_DIR = _TROPOMI_L3_ROOT / 'TROPOMI_NO2AMFtotal_01deg'
TROPOMI_01_NO2_AMF_TROP_DIR = _TROPOMI_L3_ROOT / 'TROPOMI_NO2AMFtrop_01deg'

# TM5-MP surface pressure, used when regridding TROPOMI AKs vertically
TM5_SURFACE_P_DIR = DATASETS_ROOT / 'TM5MP_Model'

# TROPOMI products regridded WITH air-mass factors, 0.05 deg over CONUS.
# These DO exist (34 daily files per species) and are the only AMF fields
# currently available anywhere. They are NOT what the vertical-regrid loaders
# read: those expect 0.15 deg (or 0.1 deg) GLOBAL fields in separate
# AMFtotal_/AMFtrop_ files, whereas these are 0.05 deg CONUS with all AMFs
# (trop/strat/total/clear/cloud) in one file per day. Using them for the NO2
# tropospheric-AK conversion would mean adapting the loader; that has
# deliberately not been done. Recorded here so the data can be found.
TROPOMI_005_WITHAMF_DIRS = {
    'NO2': DATASETS_ROOT / 'S5P_L2__NO2____HiR_2' / 'withAMF_CONUS_005',
    'HCHO': DATASETS_ROOT / 'S5P_L2__HCHO___HiR_2' / 'withAMF_CONUS_005',
}

# Species-keyed view of the 0.15 deg regridded TROPOMI L3 products,
# for scripts that loop over species.
TROPOMI_015_DIRS = {
    'NO2': TROPOMI_015_NO2_DIR,
    'HCHO': TROPOMI_015_HCHO_DIR,
    'CO': TROPOMI_015_CO_DIR,
}


# ============================================================
# 6. SURFACE OBSERVATIONS (AQS / SLAMS)
# ============================================================

AQS_DIR = DATASETS_ROOT / 'AQS'
AQS_2018_DIR = AQS_DIR / 'ForYear2018'
AQS_MATCHED_ROOT = AQS_DIR / 'MUSICA_Matched'

WRFCMAQ_LISTOS_DIR = DATASETS_ROOT / 'WRFCMAQ_LISTOS'
# CMAQ ACONC output used for the LISTOS model-vs-model comparison.
# NOTE: the delivered file names carry a per-recipient suffix, so callers must
# resolve a day's file by glob (CCTM_ACONC_*_<YYYYMMDD>*.nc) rather than by
# constructing an exact name.
CMAQ_ACONC_DIR = _env_path('MUSICA_ENV_CMAQ_ACONC_DIR',
                           WRFCMAQ_LISTOS_DIR / 'CMAQ_output' / 'ACONC_v1')
SITE_HOURLY_JJA2018_DIR = DATASETS_ROOT / 'Sitei_Hourly_JJA2018'

# CSV files mapping AQS monitor IDs to SE column indices.
# These used to live under a personal ~/Scripts/ tree, which was deleted on
# 2026-08-31; they now sit with the other processed data on the data volume.
# Written by model_evaluation/GetMatched_*_ColumnIndex.py and
# model_evaluation/func_MUSICAvsSLAMS_surf.py; read by the Match{Daily,Hourly}
# SLAMS scripts. Small (~580 KB) and regenerable, but only by re-running the
# matching against the AQS data and a spin-up h1 file.
MUSICA_COLIDX_DIR = _env_path('MUSICA_ENV_COLIDX_DIR',
                              MUSICA_PROJECT_DATA / 'AQS_colidxCSV')


def aqs_matched_dir(casename, freq='hourly'):
    """Directory of model-matched AQS files for `casename` ('hourly' or 'daily')."""
    return AQS_MATCHED_ROOT / casename / freq


# ============================================================
# 7. EMISSIONS
# ============================================================

_EMIS_ROOT = NCAR_COPIES_ROOT / 'acom' / 'MUSICA' / 'emissions'

# Biomass burning (QFED 2.6 / FINN)
BB_EMIS_NE0CONUS_DIR = _EMIS_ROOT / 'qfed2.6_finn' / 'ne0conus30x8'
BB_EMIS_NE30NP4_DIR = _EMIS_ROOT / 'qfed2.6_finn' / 'ne30np4'

# Anthropogenic (CAMS-GLOB-ANT)
CAMS_V51_ORIG_DIR = _EMIS_ROOT / 'cams' / 'CAMS-GLOB-ANTv5.1' / 'CAMS-GLOB-ANT_v5.1_orig'
CAMS_V62_ORIG_DIR = _EMIS_ROOT / 'cams' / 'CAMS-GLOB-ANT_v6.2' / 'CAMS-GLOB-ANT_v6.2_orig'
CAMS_V62_NE30NP4_DIR = _EMIS_ROOT / 'cams' / 'CAMS-GLOB-ANT_v6.2' / 'ne30np4'
CAMS_V62_NE30NP4_NONAN_DIR = _EMIS_ROOT / 'cams' / 'CAMS-GLOB-ANT_v6.2' / 'nonNANne30np4'

# NEI preprocessed to the 0.1 deg CONUS grid
NEI2017_01_DIR = CESM22_ROOT / 'CAMS_withCONUS2017NEI' / 'NEI2017_CONUS_output_01deg'
NEI2022V2_01_DIR = CESM22_ROOT / 'CAMS6.2_withCONUS2022v2NEI' / 'NEI2022v2_T1_CONUS_output_01deg'

# CAMS + NEI merged at 0.1 deg global
CAMS_NEI2017_MERGED_ROOT = CESM22_ROOT / 'CAMS_withCONUS2017NEI' / 'GLOB_Merged_CAMS_NEI_01deg'
CAMS_NEI2022_MERGED_ROOT = CESM22_ROOT / 'CAMS6.2_withCONUS2022v2NEI' / 'globCAMS_conusNEI_01deg'


def cams_nei_merged_hourly_dir(year, root=None):
    """Per-year hourly merged CAMS+NEI directory."""
    return Path(root or CAMS_NEI2022_MERGED_ROOT) / f'{year}_hourly_timeFixed'


def cams_nei_merged_byspecies_dir(year, root=None):
    """Per-year species-grouped merged CAMS+NEI directory."""
    return Path(root or CAMS_NEI2022_MERGED_ROOT) / f'{year}_GroupBySpecies'


# ============================================================
# 8. PROJECT: CONUSBGO3  —  CONUS background ozone (ne30np4)
#    Emission-zeroing experiments; Apr-Oct 2022 and 2023.
# ============================================================

# Directory names below are the existing on-disk locations; override with the
# env vars if your copy is laid out differently.
BGO3_ROOT = _env_path('MUSICA_ENV_BGO3_ROOT',
                      PROCESSED_DATA_ROOT / 'CONUSBGO3_postprocessing')

BGO3_MONITOR_INFO_DIR = BGO3_ROOT / 'MonitorInfo'
BGO3_MONITOR_LIST = _env_path('MUSICA_ENV_BGO3_MONITOR_LIST',
                              BGO3_MONITOR_INFO_DIR / 'given_monitors.csv')
BGO3_MONITOR_COLIDX = BGO3_MONITOR_INFO_DIR / 'MatchedMonitors_ne30_ColIdx.csv'

BGO3_MERGED_SURFO3_DIR = BGO3_ROOT / 'h2_surfO3_merged'
BGO3_GIVEN_MONITORS_DIR = BGO3_ROOT / 'ForGivenMonitors'
BGO3_REGRIDDED_1DEG_DIR = BGO3_ROOT / 'Regridded1deg'
BGO3_REGRIDDED_1DEG_HOURLY_DIR = BGO3_REGRIDDED_1DEG_DIR / 'hourly'

BGO3_FIGURE_DIR = FIGURES_ROOT / 'CESM_analysis' / 'BGO3'

# --- experiment matrix -------------------------------------------------
# Keyed by (scenario, year); the case name is the CESM case directory name.
BGO3_CASE_PREFIX = 'f.e22.FCnudged.ne30_ne30_mg17.BGO3.'

BGO3_CASES = {
    ('BASE', 2022): BGO3_CASE_PREFIX + 'BASEY20220401TY20230401',
    ('BASE', 2023): BGO3_CASE_PREFIX + 'BASEY20230401TY20231101',
    ('noAnthro', 2022): BGO3_CASE_PREFIX + 'noANTHROemisCONUS80kmBufferY20220401TY20221101',
    ('noAnthro', 2023): BGO3_CASE_PREFIX + 'noANTHROemisCONUS80kmBufferY20230401TY20231101',
    ('noBB', 2022): BGO3_CASE_PREFIX + 'noBBemisCONUS80kmBufferY20220401TY20221101',
    ('noBB', 2023): BGO3_CASE_PREFIX + 'noBBemisCONUS80kmBufferY20230401TY20231101',
}

BGO3_SCENARIOS = ['BASE', 'noAnthro', 'noBB']
BGO3_YEARS = [2022, 2023]

# Short label per case, e.g. 'BASE2022'; used in output file names.
BGO3_CASE_LABELS = {case: f'{scen}{yr}' for (scen, yr), case in BGO3_CASES.items()}

# Ozone-season window (MM-DD, inclusive) analysed for every case
BGO3_START_MMDD = '04-01'
BGO3_END_MMDD = '10-31'


def bgo3_unified_mda8_glob():
    """Glob pattern for the unified gridded-MDA8 deliverable.

    The provenance tag and creation date are part of the file name, so callers
    glob and take the newest match rather than assuming an exact name.
    """
    return str(BGO3_REGRIDDED_1DEG_DIR /
               'MUSICAv0_ne30_CONUS1x1_MDA8O3_BGO3_2022-2023_AprOct_*.nc')


def bgo3_merged_surfo3_glob(casename):
    """Glob pattern for the merged hourly surface-O3 file of `casename`.

    The date span is part of the file name and differs per case, so callers
    glob and take the newest match rather than assuming an exact name.
    """
    return str(BGO3_MERGED_SURFO3_DIR / f'{casename}.cam.h2.surflev.O3.*.nc')


# ============================================================
# 9. OUTPUT PROVENANCE
#    Written into NetCDF global attributes and output file names.
#    No personal name lives here; override with env vars if you want a
#    stable tag across users (e.g. an institutional one).
# ============================================================

AUTHOR_TAG = os.environ.get('MUSICA_ENV_AUTHOR_TAG', USER)
PROCESSED_BY = os.environ.get('MUSICA_ENV_PROCESSED_BY', USER)
CONTACT = os.environ.get('MUSICA_ENV_CONTACT', '')
INSTITUTION = os.environ.get(
    'MUSICA_ENV_INSTITUTION',
    'Smithsonian Astrophysical Observatory '
    '(Center for Astrophysics | Harvard & Smithsonian)')
MACHINE = os.environ.get('MUSICA_ENV_MACHINE', 'MIT Svante (svante9.mit.edu)')


# ============================================================
# 10. HELPERS
# ============================================================

def ensure_dir(path):
    """Create `path` (and parents) if needed; return it as a Path."""
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p
