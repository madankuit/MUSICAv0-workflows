#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_gridded_MDA8_3scenarios.py

Three-panel CONUS figure of the unified 1x1 gridded MDA8 O3 deliverable
(produced by Regrid_ne30_surfO3_to_1x1_conserve.py): BASE, noAnthro and noBB
shown as **absolute** seasonal-mean MDA8 over Apr-Oct 2022-2023, all on one
shared colour scale so the scenarios can be compared directly.

Companion to plot_gridded_MDA8_scenarios.py, which instead shows BASE plus the
two *differences* on independent diverging scales. Use this one to see the
absolute fields; use that one to see the emission contributions.

The colour range is set from the 2nd/98th percentiles of all three fields
together, so no scenario is clipped relative to the others.

Run in the `base` env (matplotlib + cartopy). All paths come from
config/paths.py; the deliverable is located by glob, so no file name is
hard-coded here.

MODIFICATION HISTORY:
    8 Jul 2026: VERSION 1.0
    - Initial version
    31 Aug 2026: VERSION 1.1
    - Paths moved to config/paths.py; deliverable located by glob
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
from config.paths import BGO3_FIGURE_DIR, bgo3_unified_mda8_glob, ensure_dir

_hits = sorted(glob.glob(bgo3_unified_mda8_glob()))
if not _hits:
    raise FileNotFoundError(
        "No unified MDA8 file found; run Regrid_ne30_surfO3_to_1x1_conserve.py first.")
MDA8_FILE = _hits[-1]
FIG = str(ensure_dir(BGO3_FIGURE_DIR / "regrid_trial")) + "/"
EXT = [-125, -66, 23, 50]
# ============================================================

m = xr.open_dataset(MDA8_FILE)
M = m["MDA8O3"]; lon = m.lon.values; lat = m.lat.values
names = [str(s) for s in m.scenario.values]                      # BASE, noAnthro, noBB
means = [M.isel(scenario=i).mean("time") for i in range(len(names))]   # 2-yr Apr-Oct mean

# One shared scale across all three scenarios
allv = np.concatenate([mn.values.ravel() for mn in means])
vmin = float(np.floor(np.nanpercentile(allv, 2)))
vmax = float(np.ceil(np.nanpercentile(allv, 98)))
print("source:", MDA8_FILE)
print("vmin/vmax:", vmin, vmax)
for n, mn in zip(names, means):
    print(f"  {n:9s} 2yr-mean MDA8 ppb: mean={float(mn.mean()):.1f} "
          f"min={float(mn.min()):.1f} max={float(mn.max()):.1f}")

fig, axs = plt.subplots(1, 3, figsize=(16, 4.9),
                        subplot_kw={"projection": ccrs.PlateCarree()})
for ax, mn, n in zip(axs, means, names):
    pc = ax.pcolormesh(lon, lat, mn, cmap="YlOrRd", vmin=vmin, vmax=vmax,
                       shading="auto", transform=ccrs.PlateCarree())
    ax.add_feature(cfeat.COASTLINE, lw=.5)
    ax.add_feature(cfeat.STATES, lw=.3, edgecolor="gray")
    ax.add_feature(cfeat.BORDERS, lw=.4)
    ax.set_extent(EXT, ccrs.PlateCarree())
    ax.set_title(n, fontsize=13, fontweight="bold")

cb = fig.colorbar(pc, ax=axs, orientation="horizontal", shrink=0.45,
                  pad=0.06, aspect=45, extend="both")
cb.set_label("seasonal-mean MDA8 O3 (ppb)", fontsize=11)
fig.suptitle("MUSICAv0 ne30 (1x1, mass-conservative) — seasonal-mean MDA8 O3, "
             "Apr-Oct 2022-2023", fontsize=14, y=0.98)

out = FIG + "MDA8_3scenarios_seasonmean_2022-2023.png"
fig.savefig(out, dpi=140, bbox_inches="tight")
print("saved", out)
