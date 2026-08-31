#!/usr/bin/env python
"""Generate ne30np4 -> 1x1 FV CONSERVATIVE ESMF weights (esmpy 8.7),
mirroring the Regridding_ESMF_v1 call_ESMF setup. Weight-file only + verify."""
import esmpy, numpy as np, xarray as xr, datetime, os
import sys, pathlib

# Paths come from config/paths.py - the single source of truth.
_ROOT = next(p for p in pathlib.Path(__file__).resolve().parents
             if (p / 'config' / 'paths.py').exists())
sys.path.insert(0, str(_ROOT))
import config  # noqa: F401
from config.paths import SCRIP_NE30NP4, FV_GRIDINFO_1X1, GRIDS_EXTERNAL_DIR, ensure_dir

SRC = str(SCRIP_NE30NP4)      # SE source SCRIP (corners), ships with the repo
DST = str(FV_GRIDINFO_1X1)    # 1x1 global FV grid info (CF bnds), ships with the repo
YMD = datetime.datetime.now().strftime("%Y%m%d")
# The weight file is large and derived, so it is written outside the repo.
WGT = str(ensure_dir(GRIDS_EXTERNAL_DIR) / f"ESMFmap_ne30np4_TO_1x1_conserve_c{YMD}.nc")

print("esmpy", esmpy.__version__)
# --- source: SE mesh from SCRIP, field on elements ---
srcmesh  = esmpy.Mesh(filename=SRC, filetype=esmpy.FileFormat.SCRIP)
srcfield = esmpy.Field(srcmesh, name="src", meshloc=esmpy.MeshLoc.ELEMENT)
# --- destination: FV grid from GRIDSPEC, field on centers ---
dstgrid  = esmpy.Grid(filename=DST, filetype=esmpy.FileFormat.GRIDSPEC,
                      add_corner_stagger=True)
dstfield = esmpy.Field(dstgrid, name="dst", staggerloc=esmpy.StaggerLoc.CENTER)

print("generating CONSERVE weights ->", WGT)
t0 = datetime.datetime.now()
rg = esmpy.Regrid(srcfield, dstfield, filename=WGT,
                  regrid_method=esmpy.RegridMethod.CONSERVE,
                  unmapped_action=esmpy.UnmappedAction.IGNORE,
                  ignore_degenerate=True)
print("weight-gen time:", datetime.datetime.now() - t0)

# --- verify conservative + orientation ---
ds = xr.open_dataset(WGT)
col, row, S = ds["col"].values, ds["row"].values, ds["S"].values
print("n_s:", len(S))
print("col(SOURCE) max:", col.max(), "(ne30np4=48602)")
print("row(DEST)   max:", row.max(), "(1x1 =181*361=", 181*361, ")")
uniq, idx = np.unique(np.sort(row), return_index=True)
sums = np.add.reduceat(S[np.argsort(row)], idx)
print(f"dest weight-sum over MAPPED cells: min={sums.min():.4f} "
      f"mean={sums.mean():.4f} max={sums.max():.4f}  (mapped {uniq.size}/{181*361})")
print("SAVED:", WGT)
