#!/usr/bin/env python
"""Generate ne30np4 -> 1x1 FV CONSERVATIVE ESMF weights (esmpy 8.7),
mirroring the Regridding_ESMF_MTv1 call_ESMF setup. Weight-file only + verify."""
import esmpy, numpy as np, xarray as xr, datetime, os

G = "/net/fs09/d0/taoma528/CESM22/grids/"
SRC = G + "ne30np4_091226_pentagons.nc"      # SE source SCRIP (corners)
DST = G + "FV1x1grid_info_c20241105.nc"       # 1x1 global FV grid info (CF bnds)
YMD = datetime.datetime.now().strftime("%Y%m%d")
WGT = G + f"ESMFmap_ne30np4_TO_1x1_conserve_c{YMD}.nc"

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
