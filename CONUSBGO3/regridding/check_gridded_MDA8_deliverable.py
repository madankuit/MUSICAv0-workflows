#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
check_gridded_MDA8_deliverable.py

Quality-control the unified 1x1 gridded MDA8 O3 deliverable before it is shared.

Two things happen here:

1. **Coverage report** — for every (scenario, year) combination, how many daily
   dates are present and what fraction of cells are NaN. This is the check that
   surfaced the known local-time edge case: Oct 31 in the western time zones is
   NaN for BASE-2023 and noBB-2022, because those merged inputs stop at
   Nov 1 00 UTC and the last local day therefore has fewer than 13 of the 17
   valid 8-h windows, so it is correctly dropped (~0.2 % of values). A specific
   west-coast cell is printed so that case stays visible.

2. **Per-year deliverable maps** — BASE seasonal mean, plus BASE-noAnthro and
   BASE-noBB, for a single year (default 2022). These are the "does it look
   right" figures, distinct from the two-year publication figures produced by
   plot_gridded_MDA8_scenarios.py and plot_gridded_MDA8_3scenarios.py.

Run in the `base` env (matplotlib + cartopy). All paths come from
config/paths.py; the deliverable is located by glob, so no file name is
hard-coded here.

MODIFICATION HISTORY:
    8 Jul 2026: VERSION 1.0
    - Initial version
    31 Aug 2026: VERSION 1.1
    - Paths moved to config/paths.py; deliverable located by glob; the checked
      year is a module constant instead of being hard-coded throughout
"""
import glob, sys, pathlib
import numpy as np, xarray as xr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import cartopy.crs as ccrs, cartopy.feature as cfeat

# ============================================================
# CONFIGURATION - every path comes from config/paths.py.
# ============================================================
_ROOT = next(p for p in pathlib.Path(__file__).resolve().parents
             if (p / 'config' / 'paths.py').exists())
sys.path.insert(0, str(_ROOT))
import config  # noqa: F401
from config.paths import (BGO3_FIGURE_DIR, BGO3_YEARS,
                          bgo3_unified_mda8_glob, ensure_dir)

_hits = sorted(glob.glob(bgo3_unified_mda8_glob()))
if not _hits:
    raise FileNotFoundError(
        "No unified MDA8 file found; run Regrid_ne30_surfO3_to_1x1_conserve.py first.")
MDA8_FILE = _hits[-1]
FIG = str(ensure_dir(BGO3_FIGURE_DIR / "regrid_trial")) + "/"

CHECK_YEAR = 2022          # year mapped in the per-year figures
EXT = [-125, -66, 23, 50]
# a west-coast (PDT) cell, to keep the Oct-31 local-time edge case visible
EDGE_LAT, EDGE_LON, EDGE_DATE = 34, -118, "2023-10-31"
# ============================================================

m = xr.open_dataset(MDA8_FILE)
M = m["MDA8O3"]
print("source:", MDA8_FILE, "\n")

# ---------- 1. coverage ----------
print("coverage (days present, and % of cells NaN):")
for si, sc in enumerate(m.scenario.values):
    for yr in BGO3_YEARS:
        sub = M.isel(scenario=si).sel(time=str(yr))
        print(f"  {str(sc):9s} {yr}: days={sub.sizes['time']:3d} "
              f"nan%={100 * float(np.isnan(sub).mean()):.2f}")

w = M.isel(scenario=0).sel(time=EDGE_DATE).sel(
    lat=EDGE_LAT, lon=EDGE_LON, method="nearest")
print(f"\nBASE {EDGE_DATE} @ ({EDGE_LAT},{EDGE_LON}) west-coast cell MDA8 = "
      f"{float(w):.1f} ppb   (NaN here is the expected local-time edge case)\n")

# ---------- 2. per-year maps ----------
sm = {str(sc): M.isel(scenario=si).sel(time=str(CHECK_YEAR)).mean("time")
      for si, sc in enumerate(m.scenario.values)}
lon = m.lon.values; lat = m.lat.values
proj = dict(projection=ccrs.PlateCarree())


def basemap(ax):
    ax.add_feature(cfeat.COASTLINE, lw=.5)
    ax.add_feature(cfeat.STATES, lw=.3, edgecolor="gray")
    ax.add_feature(cfeat.BORDERS, lw=.4)
    ax.set_extent(EXT, ccrs.PlateCarree())


fig, ax = plt.subplots(figsize=(8, 5), subplot_kw=proj)
pc = ax.pcolormesh(lon, lat, sm["BASE"], cmap="YlOrRd", vmin=30, vmax=65,
                   shading="auto", transform=ccrs.PlateCarree())
basemap(ax)
plt.colorbar(pc, ax=ax, shrink=.8, label="MDA8 O3 (ppb)")
ax.set_title(f"BASE seasonal-mean MDA8 O3  Apr-Oct {CHECK_YEAR}  (1x1 conservative)")
plt.tight_layout()
plt.savefig(FIG + f"DELIV_BASE_MDA8_seasonmean_{CHECK_YEAR}.png", dpi=130)
plt.close()

for k in ["noAnthro", "noBB"]:
    d = sm["BASE"] - sm[k]
    vmax = float(np.nanpercentile(np.abs(d), 99))
    fig, ax = plt.subplots(figsize=(8, 5), subplot_kw=proj)
    pc = ax.pcolormesh(lon, lat, d, cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                       shading="auto", transform=ccrs.PlateCarree())
    basemap(ax)
    plt.colorbar(pc, ax=ax, shrink=.8, label="delta MDA8 O3 (ppb)")
    ax.set_title(f"BASE - {k}  seasonal-mean MDA8 O3  Apr-Oct {CHECK_YEAR}")
    plt.tight_layout()
    plt.savefig(FIG + f"DELIV_BASE_minus_{k}_MDA8_seasonmean_{CHECK_YEAR}.png", dpi=130)
    plt.close()
    print(f"BASE-{k} seasonmean MDA8: max +{np.nanmax(d):.1f} / "
          f"min {np.nanmin(d):.1f} ppb")

print("\nfigures ->", FIG)
