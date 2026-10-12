#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
diag_grid_mda8_mismatch.py

Follow-up to audit_deliverables.py part D: scan many cells of one unified MDA8 file,
recompute MDA8 (audit's numpy code) from the matching hourly gridded file, and for
every cell that differs report its location, the UTC offset the audit used, which
offset reproduces the file, and on which days the values differ.

    python diag_grid_mda8_mismatch.py ne0CONUSne30x8 0.15 noBBCONUS 2024 [ncells]

    12 Oct 2026: VERSION 1.0
"""
import sys, pathlib, glob
import numpy as np, pandas as pd, xarray as xr

HERE = pathlib.Path(__file__).resolve().parent
src = (HERE.parent / 'audit_deliverables.py').read_text()
exec(src[:src.index('# ======================================================================= A')])  # helpers only

EXP, TGT, SC, YR = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
N = int(sys.argv[5]) if len(sys.argv) > 5 else 400
E = BGO3_EXPERIMENTS[EXP]; tag = E['targets'][TGT][2]; lab = f"{E['label_prefix']}{SC}{YR}"
U = xr.open_dataset(sorted(glob.glob(str(E['regrid_dir'] / (E['unified_stem'].format(tag=tag) + '_c*.nc'))))[-1]).MDA8O3
H = xr.open_dataset(sorted(glob.glob(str(E['hourly_dir'] / f'{tag}_UTChourlySurfO3.{lab}.*.nc')))[-1]).O3
rng = np.random.default_rng(1); bad = 0
for _ in range(N):
    iy = int(rng.integers(0, H.sizes['lat'])); ix = int(rng.integers(0, H.sizes['lon']))
    s = H.isel(lat=iy, lon=ix).to_series(); la, lo = float(H.lat[iy]), float(H.lon[ix])
    off = summer_offset(la, lo)
    th = U.sel(scenario=SC, time=slice(f'{YR}-04-01', f'{YR}-10-30')).isel(lat=iy, lon=ix).values
    mine = mda8_numpy(s.index, s.values, off, YR).loc[:f'{YR}-10-30']
    dd = np.abs(mine.values - th)
    if np.nanmax(dd) > 1e-3:
        bad += 1
        per = {o: round(float(np.nanmax(np.abs(mda8_numpy(s.index, s.values, o, YR).loc[:f'{YR}-10-30'].values - th))), 4)
               for o in range(-9, -3)}
        days = list(mine.index[dd > 1e-3].strftime('%m-%d'))
        print(f"cell lat {la:.3f} lon {lo:.3f}: audit offset {off}, max|d| {np.nanmax(dd):.2f} ppb on {len(days)} days "
              f"{days[:6]} | max|d| by offset {per}")
print(f"{lab} {TGT}: {N} cells checked, {bad} differ")
