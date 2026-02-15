"""
Script: ne30CONUSne30x8_scale_emissions_by_Canadaregion_c20260212.py

Description
-----------
This standalone Python script scales gridded emissions within a specified
geographic region by a user-defined multiplicative factor on the ne30CONUSne30x8 grid.

The script:
1. Reads gridded emissions data (NetCDF format).
2. Applies a spatial mask defining the target region.
3. Multiplies emissions within the masked region by a specified scaling factor.
4. Writes the modified emissions to a new output file, preserving metadata.

Inputs
------
- BBEmissions_diri: Original emissions directory (NetCDF files); Scaled files are saved to a subdirectory under the original input file directory.
- varname: Target variable name (e.g., "NO", "CO")
- regionMask: Region mask (boolean array or shapefile-derived mask); should be defined already
- scalefactor: Scaling factor (float; e.g., 0.7 for 30% reduction, 1.3 for 30% increase)

Output
------
- NetCDF files containing scaled emissions

MODIFICATION HISTORY:
    Madankui Tao, Feb, 12, 2026: VERSION 1
    - Initial version

"""

# ============================================================
# USAGE INSTRUCTIONS
# ------------------------------------------------------------
# 1) Set BBEmissions_diri at the top of this script to the
#    directory containing your original emissions files.
#
# 2) Adjust scalefactor and regionMask as needed.
#
# 3) Uncomment the execution code in the end and run the script from the terminal:
#
#       python scale_emissions_by_canada_region.py
#
# The scaled files will be written to:
#       {BBEmissions_diri}/ne0conus30x8_<region>Masked_<XXpct>/
#
# Example:
#       python scale_emissions_by_canada_region.py
# ============================================================

#================================================================================================
# inputs to change
varname = 'CO' # BB emissions variable if not to process for all ('all')

regionMask = "mask_Quebec"
scalefactor = 0.7   # e.g. 0.7=70pct, 1.3=130pct

# # Directory containing biomass-burning emissions files
# BBEmissions_diri = "path/to/biomass_burning_emissions/"
# # Path to precomputed regional mask file for ne0CONUSne30x8 grid
# MaskFile = "path/to/mask_CanadaProvinces_ne0CONUS_ne30x8.nc"

#================================================================================================
#### Module import ###
import os
import re
import glob
import fnmatch
from pathlib import Path

import pandas as pd
import geopandas as gpd
from shapely import wkt
from shapely.geometry import Point, Polygon

import xarray as xr
import numpy as np

#================================================================================================
### Get output directory
out_diri = Path(Out_BBEmissions_diri)
# Create directory if it does not exist
out_diri.mkdir(parents=True, exist_ok=True)
print(f"Output directory ready: {out_diri}")

# Find all bb emissions files
pattern = os.path.join(
    BBEmissions_diri,
    "qfed.emis_*_hires_mol_2020-2024_c20250115.nc"
)
files = glob.glob(pattern)
regex = re.compile(
    r"qfed.emis_(?P<spc>.+?)_hires_mol_2020-2024_c20250115\.nc"
)
spc_list = [
    regex.search(os.path.basename(f)).group("spc")
    for f in files
    if regex.search(os.path.basename(f))
]

# Get the file(s) to modify
if varname=='all':
    scalespc_list = spc_list
    scalefiles = files
else:
    # for a single variable
    scalespc_list = [varname]
    scalefiles = f'{BBEmissions_diri}qfed.emis_{varname}_hires_mol_2020-2024_c20250115.nc'
    
#================================================================================================
# Grid and Mask Information
mask_ds = xr.open_dataset(MaskFile)
# for the selected region mask
mask_da_loaded = mask_ds[regionMask].astype(bool)  # dims: (grid_size,) or (ncol,

# Since the emissions files use dimension name "ncol" and mask uses "grid_size",
# align by renaming mask dim to match emissions.
if "grid_size" in mask_da_loaded.dims and "ncol" not in mask_da_loaded.dims:
    mask_da_loaded = mask_da_loaded.rename({"grid_size": "ncol"})

# Coordinate-like vars to clean attributes (and fill NaN only if numeric)
coord_like_vars = {'time', 'ncol', 'lat', 'lon', 'area', 'date', 'altitude', 'rrfac'}

# where to put the output files based on the input
# Convert scaling factor to percentage label
pct_value = int(round(scalefactor * 100))   # 0.7 → 70
pct_str = f"{pct_value}pct"

BASE_OUT = Path(BBEmissions_diri)   # change locally if needed
# The scaled files will be placed in a subdirectory from the original file
Out_BBEmissions_diri = BASE_OUT / f"ne0conus30x8_{regionMask}Masked_{pct_str}"
Out_BBEmissions_diri.mkdir(parents=True, exist_ok=True)

#================================================================================================
### Key function
def scale_bb_emissions_by_mask(
    spc,
    in_diri,
    out_diri,
    mask_da,
    scalefactor,
    regionMask,
    coord_like_vars=None,
):
    """
    Scale wildfire (QFED) emissions for a given species within a specified region mask.

    Parameters
    ----------
    spc : str
        Species name (e.g., "NO", "CO").
    in_diri : str
        Directory containing original emission files.
    out_diri : str
        Directory to write scaled emission files.
    mask_da : xr.DataArray (bool)
        Boolean mask on the same 'ncol' grid (True = region to scale).
    scalefactor : float
        Multiplicative factor applied inside the mask.
    regionMask : str
        Name of the mask variable (for metadata documentation).
    coord_like_vars : iterable, optional
        Variables to exclude from scaling (coordinates, metadata fields).

    Notes
    -----
    - Scaling is applied only to variables containing the 'ncol' dimension.
    - Outside the mask, values remain unchanged.
    - Border classification follows cell-center logic.
    """

    import os
    import numpy as np
    import xarray as xr

    if coord_like_vars is None:
        coord_like_vars = set()

    try:
        # --------------------------------------------------
        # Open input emissions file
        # --------------------------------------------------
        spc_fileINpath = (
            f"{in_diri}qfed.emis_{spc}_hires_mol_2020-2024_c20250115.nc"
        )
        spci_ds = xr.open_dataset(spc_fileINpath)

        # --------------------------------------------------
        # Optional: adjust longitudes (0–360 → -180–180)
        # Diagnostic only; does not affect scaling
        # --------------------------------------------------
        if "lon" in spci_ds:
            lon = spci_ds["lon"].values
            adj_lon = np.where(lon >= 180, lon - 360, lon)
            spci_ds["Adjustedlons"] = xr.DataArray(
                adj_lon,
                dims=("ncol",),
                coords={"ncol": spci_ds["ncol"].values},
            )

        # --------------------------------------------------
        # Apply regional scaling
        # inside mask:  value * scalefactor
        # outside mask: unchanged
        # --------------------------------------------------
        for var in spci_ds.data_vars:

            if var in coord_like_vars:
                continue

            da = spci_ds[var]

            # Only scale variables defined on the unstructured grid
            if "ncol" not in da.dims:
                continue

            spci_ds[var] = da.where(~mask_da, other=da * scalefactor)

            # Record scaling metadata
            spci_ds[var].attrs.update({
                "region_scaling_applied": "True",
                "region_mask_name": regionMask,
                "region_scale_factor": float(scalefactor),
            })

        # Remove temporary lon variable if created
        if "Adjustedlons" in spci_ds:
            spci_ds = spci_ds.drop_vars("Adjustedlons")

        # --------------------------------------------------
        # Write output file
        # --------------------------------------------------
        os.makedirs(out_diri, exist_ok=True)

        spc_fileOUTpath = (
            f"{out_diri}qfed.emis_{spc}_hires_mol_2020-2024_c20250115.nc"
        )

        spci_ds.to_netcdf(spc_fileOUTpath)

        print(f"Saved scaled emissions: {spc_fileOUTpath}")
        print("*" * 72)

    except Exception as e:
        print(f"Error encountered for {spc}: {e}")

# #================================================================================================
# ### Execution Step
# # ----------------------------
# # Scale emissions within mask
# # ----------------------------
# for spc in scalespc_list:
#     scale_bb_emissions_by_mask(
#         spc,
#         BBEmissions_diri,
#         Out_BBEmissions_diri,
#         mask_da_loaded,
#         scalefactor,
#         regionMask,
#         coord_like_vars,
#     )