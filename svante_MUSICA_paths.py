"""
svante_MUSICA_paths.py

Reference file for all MUSICA-related data paths on Svante.
Copy the relevant variables into your script's USER CONFIGURATION block,
or import with:  from svante_MUSICA_paths import *

All paths are on the Svante file system (/net/fs09/d0/taoma528/).
Local script paths are on the Svante home directory (/home/taoma528/).

MODIFICATION HISTORY:
    Madankui Tao, Mar, 2026: VERSION 1.0
    - Initial version, compiled from MUSICAv0-workflows scripts
"""

# ============================================================
# ROOT DIRECTORIES
# ============================================================

svante_archive    = '/net/fs09/d0/taoma528/CESM22/archive/'
CESM22_root       = '/net/fs09/d0/taoma528/CESM22/'
Datasets_root     = '/net/fs09/d0/taoma528/Datasets/'
ProcessedData_root = '/net/fs09/d0/taoma528/ProcessedData/'
ncar_copies_root  = '/net/fs09/d0/taoma528/ncar_copies/'

# ============================================================
# GRID FILES (SCRIP & ESMF regridding weights)
# ============================================================

# SCRIP files — grid geometry descriptions
SCRIP_ne0CONUSne30x8 = '/net/fs09/d0/taoma528/CESM22/grids/ne0CONUS_ne30x8_np4_SCRIP.nc'
SCRIP_ne30np4        = '/net/fs09/d0/taoma528/CESM22/grids/ne30np4_091226_pentagons.nc'
SCRIP_ne30np4_new    = '/net/fs09/d0/taoma528/CESM22/grids/ne30np4_grid_c20241105.nc'

# FV target grid info files (needed as ESMF destination grid)
ESMFmap_015grid_file = '/net/fs09/d0/taoma528/CESM22/grids/FV_gridinfo_0.15_c20231204.nc'
ESMFmap_01grid_file  = '/net/fs09/d0/taoma528/CESM22/grids/FV_gridinfo_CAMS_c20210219.nc'
ESMFmap_1x1grid_file = '/net/fs09/d0/taoma528/CESM22/grids/FV1x1grid_info_c20241105.nc'

# ESMF conservative regridding weight files
Regridding_015weights_ne0CONUS  = '/net/fs09/d0/taoma528/CESM22/grids/ne0CONUSne30x8_ESMFmap_0.15x0.15_cubit_conserve_cams_c20231204.nc'
Regridding_01weights_ne0CONUS   = '/net/fs09/d0/taoma528/CESM22/grids/ESMFmap_0.1x0.1_ne0CONUSne30x8_cubit_conserve_cams.nc'
Regridding_09x125weights_ne30   = '/net/fs09/d0/taoma528/CESM22/grids/ESMFmap_0.9x1.25_ne30np4_cubit_conserve_cams.nc'
Regridding_1x1weights_ne30      = '/net/fs09/d0/taoma528/CESM22/grids/ESMFmap_0.9x1.25_ne30np4_cubit_conserve_cams.nc'

# ============================================================
# CESM/MUSICA OUTPUT — RAW HISTORY FILES
# ============================================================

# Raw h1/h2 history files live under:
#   svante_archive + casename + '/atm/hist/'
# Example:
#   casename = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1base01'
#   RunPath  = svante_archive + casename + '/atm/hist/'

# ============================================================
# CESM/MUSICA OUTPUT — PROCESSED / MERGED
# ============================================================

# Surface-layer merged h2 files (July 2018, ne0CONUSne30x8)
processed_output_diri         = '/net/fs09/d0/taoma528/CESM22/processed_output/'
July2018_surfh2_merged_diri   = '/net/fs09/d0/taoma528/CESM22/processed_output/July2018_surfh2_merged/'

# Vertical column density (VCD) outputs — on native ncol grid
# (computed by postprocessing/CalcTropVCD_* and CalcTotalVCD_*)
MUSICATropVCDSave_diri = '/net/fs09/d0/taoma528/CESM22/Calculated_MUSICA_VCD/TroposphericVCD/{casename}/'
MUSICAVCDSave_diri     = '/net/fs09/d0/taoma528/CESM22/Calculated_MUSICA_VCD/TotalVCD/{casename}/'

# Met file used with VCD v2 scripts (external met)
MetPath = '/net/fs09/d0/taoma528/CESM22/Calculated_MUSICA_VCD/TroposphericVCD/met_f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1base01.nc'

# ============================================================
# CESM/MUSICA OUTPUT — REGRIDDED TO 0.15° LAT-LON
# ============================================================

# Root for 2018 July TROPOMI-time regridded outputs
Regridded_MUSICA_2018_root = '/net/fs09/d0/taoma528/CESM22/Regridded_MUSICA_Output/2018_1330LT_TROPOMIcomp/MassConserve_latlon015_MUSICAoutput/'

# Per-casename subdirectory structure (append casename):
#   .../h2/                        — full hourly regridded files
#   .../h2_ByVar/                  — split by variable
#   .../h2_VertiVarsByVarDate/     — split by variable and date (vertical vars)
#   .../TropVCD_approx1330LT/      — tropospheric VCD with TROPOMI AK applied
#   .../TotalVCD_approx1330LT/     — total VCD (CO) with TROPOMI AK applied
#
# Example:
#   regridByVar_diri      = Regridded_MUSICA_2018_root + casename + '/h2_ByVar/'
#   MUSICATropVCDSave_diri = Regridded_MUSICA_2018_root + casename + '/TropVCD_approx1330LT/'
#   MUSICATotalVCDSave_diri = Regridded_MUSICA_2018_root + casename + '/TotalVCD_approx1330LT/'

# ============================================================
# TROPOMI SATELLITE DATA
# ============================================================

# L2 offline reprocessed files — native granule format (2018)
TROPOMI_NO2_L2_diri  = '/net/fs09/d0/taoma528/Datasets/S5P_L2__NO2____HiR_2/2018/'
TROPOMI_HCHO_L2_diri = '/net/fs09/d0/taoma528/Datasets/S5P_L2__HCHO___HiR_2/2018/'
TROPOMI_CO_L2_diri   = '/net/fs09/d0/taoma528/Datasets/S5P_L2__CO_____HiR_2/2018/'

# TROPOMI regridded to 0.15° (matching MUSICA 0.15° grid, 2018)
TROPOMI_NO2_015_diri  = '/net/fs09/d0/taoma528/Datasets/S5P_L2__NO2____HiR_2/regrid_2018_MUSICA015/'
TROPOMI_HCHO_015_diri = '/net/fs09/d0/taoma528/Datasets/S5P_L2__HCHO___HiR_2/regrid_2018_MUSICA015/'
TROPOMI_CO_015_diri   = '/net/fs09/d0/taoma528/Datasets/S5P_L2__CO_____HiR_2/regrid_2018_MUSICA015/'

# TROPOMI averaging kernels regridded to 0.1° global
TROPOMI_NO2AK_01_diri  = '/net/fs09/d0/taoma528/Datasets/TROPOMI_RegridL2_ToL3/Global/TROPOMI_NO2AK_01deg/'
TROPOMI_HCHOAK_01_diri = '/net/fs09/d0/taoma528/Datasets/TROPOMI_RegridL2_ToL3/Global/TROPOMI_HCHOAK_01deg/'
TROPOMI_NO2_01_diri    = '/net/fs09/d0/taoma528/Datasets/TROPOMI_RegridL2_ToL3/Global/TROPOMI_NO2_01deg/'

# TM5-MP surface pressure (used for vertical regridding of TROPOMI AKs)
TM5SurfaceP_diri = '/net/fs09/d0/taoma528/Datasets/TM5MP_Model/'

# ============================================================
# AQS / SLAMS SURFACE OBSERVATIONS
# ============================================================

AQS_diri      = '/net/fs09/d0/taoma528/Datasets/AQS/'
AQS2018_diri  = '/net/fs09/d0/taoma528/Datasets/AQS/ForYear2018/'

# Model-matched AQS files (output of MatchHourly/MatchDaily scripts)
#   Matched_diri = f'/net/fs09/d0/taoma528/Datasets/AQS/MUSICA_Matched/{casename}/hourly/'
#   Matched_diri = f'/net/fs09/d0/taoma528/Datasets/AQS/MUSICA_Matched/{casename}/daily/'

# ============================================================
# EMISSIONS
# ============================================================

# Biomass burning — QFED2.6 / FINN (ne0CONUSne30x8 grid)
BBEmissions_ne0CONUS_diri = '/net/fs09/d0/taoma528/ncar_copies/acom/MUSICA/emissions/qfed2.6_finn/ne0conus30x8/'

# Anthropogenic — CAMS v5.1 original files (used for 2018 simulations)
CAMS_v51_orig_diri = '/net/fs09/d0/taoma528/ncar_copies/acom/MUSICA/emissions/cams/CAMS-GLOB-ANTv5.1/CAMS-GLOB-ANT_v5.1_orig/'

# Anthropogenic — CAMS v6.2 original files (most current, used for 2022-2023)
CAMS_v62_orig_diri = '/net/fs09/d0/taoma528/ncar_copies/acom/MUSICA/emissions/cams/CAMS-GLOB-ANT_v6.2/CAMS-GLOB-ANT_v6.2_orig/'

# NEI 2017 preprocessed to 0.1° CONUS grid (for 2018 simulations)
NEI2017_01_diri = '/net/fs09/d0/taoma528/CESM22/CAMS_withCONUS2017NEI/NEI2017_CONUS_output_01deg/'

# NEI 2022v2 preprocessed to 0.1° CONUS grid (for 2022-2023 simulations)
NEI2022v2_01_diri = '/net/fs09/d0/taoma528/CESM22/CAMS6.2_withCONUS2022v2NEI/NEI2022v2_T1_CONUS_output_01deg/'

# CAMS + NEI merged at 0.1° global (2018 simulations)
CAMS_NEI2017_merged_root = '/net/fs09/d0/taoma528/CESM22/CAMS_withCONUS2017NEI/GLOB_Merged_CAMS_NEI_01deg/'
# Per-year hourly: CAMS_NEI2017_merged_root + '{year}_hourly/'
# Grouped by species: CAMS_NEI2017_merged_root + '{year}_GroupBySpecies/'

# CAMS v6.2 + NEI 2022v2 merged at 0.1° global (2022-2023 simulations)
CAMS_NEI2022_merged_root = '/net/fs09/d0/taoma528/CESM22/CAMS6.2_withCONUS2022v2NEI/globCAMS_conusNEI_01deg/'
# Per-year hourly time-fixed: CAMS_NEI2022_merged_root + '{year}_hourly_timeFixed/'
# Grouped by species: CAMS_NEI2022_merged_root + '{year}_GroupBySpecies/'

# ============================================================
# PROJECT-SPECIFIC OUTPUT DIRECTORIES
# ============================================================

# Dan Jaffe / background O3 project
DanJaffe_output_diri = '/net/fs09/d0/taoma528/ProcessedData/DanJaffeMUSICAPostprocessing/'

# ============================================================
# LOCAL SCRIPT / INDEX FILES (on Svante home)
# ============================================================

# CSV files mapping AQS monitor IDs to ne0CONUSne30x8 ncol indices
MUSICA_index_diri = '/home/taoma528/Scripts/CESM_analysis/ACP_MUSICANEI_scripts/colidxCSV/'
