#!/usr/bin/env python3
"""
trial_regrid_linear.py  (ARCHIVE - rejected alternative, 8 Jul 2026)

Single-day trial of a **linear (Delaunay) interpolation** regrid of ne30np4
surface O3 to a 1 degree CONUS box, for BASE/noAnthro/noBB on 2022-07-20, with
the same monitor cross-check and figures as the conservative trial.

This approach was **not** adopted. Linear interpolation between SE column
centres does not conserve mass, which matters for a deliverable meant to be
integrated and compared against gridded products; the conservative ESMF remap
in apply_conserve_trial.py was used instead and became
../Regrid_ne30_surfO3_to_1x1_conserve.py.

Kept only as the record of that comparison. Do not use it to produce results.

MODIFICATION HISTORY:
    8 Jul 2026: VERSION 1.0
    31 Aug 2026: VERSION 1.1 - paths moved to config/paths.py
"""
import time as _t
import numpy as np, pandas as pd, xarray as xr
from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================================================
# CONFIGURATION - every path comes from config/paths.py, the single
# source of truth. Override cluster locations with the MUSICA_ENV_*
# environment variables documented there. Do not hard-code paths here.
# ============================================================
import sys, pathlib, glob
_ROOT = next(p for p in pathlib.Path(__file__).resolve().parents
             if (p / 'config' / 'paths.py').exists())
sys.path.insert(0, str(_ROOT))
import config  # noqa: F401
from config.paths import (
    BGO3_CASES, BGO3_MONITOR_COLIDX, BGO3_FIGURE_DIR,
    bgo3_merged_surfo3_glob, case_hist_dir, ensure_dir,
)

DAY = "2022-07-20"          # single trial day
FIG = str(ensure_dir(BGO3_FIGURE_DIR / "regrid_trial")) + "/"


def _merged(scenario, year=2022):
    """Newest merged hourly surface-O3 file for one scenario/year."""
    hits = sorted(glob.glob(bgo3_merged_surfo3_glob(BGO3_CASES[(scenario, year)])))
    if not hits:
        raise FileNotFoundError(
            f"No merged surface-O3 file for {scenario} {year}; "
            f"run Merge_h2files_hourlysurfO3.py first.")
    return hits[-1]


merged = {k: _merged(k) for k in ("BASE", "noAnthro", "noBB")}

# Any h2 file of the BASE 2022 case provides the ne30 column coordinates.
_hist = case_hist_dir(BGO3_CASES[("BASE", 2022)])
_ex = sorted(glob.glob(str(_hist / "*.cam.h2.*.nc")))
if not _ex:
    raise FileNotFoundError(f"No h2 history files under {_hist}")
h2ex = _ex[0]
# ============================================================

# target grid (matches rewrite_ne30_1latlongrid_forCONUS_surflayer)
res=1.0
lon2d=np.arange(210,310+res,res); lat2d=np.arange(0,70+res,res)
X,Y=np.meshgrid(lon2d,lat2d)

# ne30 column lat/lon (same ncol order for every case)
llat=xr.open_dataset(h2ex)['lat'].values; llon=xr.open_dataset(h2ex)['lon'].values
t0=_t.time(); tri=Delaunay(np.column_stack([llon,llat])); print(f"Delaunay: {_t.time()-t0:.1f}s")

def regrid_day(path):
    ds=xr.open_dataset(path)
    sub=ds.sel(time=slice(DAY+"T00", DAY+"T23"))['O3']
    hrs=sub.sizes['time']; vals=sub.values*1e9  # ppb
    out=np.full((hrs,lat2d.size,lon2d.size),np.nan)
    t0=_t.time()
    for i in range(hrs):
        out[i]=LinearNDInterpolator(tri,vals[i])(X,Y)
    return out, sub.time.values, (_t.time()-t0)/hrs*1000

reg={}; nhr=None
for k,p in merged.items():
    o,tt,mspp=regrid_day(p); reg[k]=o; nhr=o.shape[0]
    print(f"{k:9s}: {o.shape}  {mspp:.0f} ms/step  dailymean ppb "
          f"min={np.nanmin(o.mean(0)):.1f} max={np.nanmax(o.mean(0)):.1f}")

# cross-check BASE gridded vs raw point column (hour 18 UTC)
md=pd.read_csv(BGO3_MONITOR_COLIDX)
md=md[md['MUSICA0_colIndex']!='Find None'].dropna(subset=['MUSICA0_colIndex'])
md['MUSICA0_colIndex']=md['MUSICA0_colIndex'].astype(int)
raw=xr.open_dataset(merged['BASE']).sel(time=slice(DAY+"T00",DAY+"T23"))['O3'].values*1e9
gi=xr.DataArray(reg['BASE'][18],dims=['lat','lon'],coords={'lat':lat2d,'lon':lon2d})
print("\ncross-check hour18  raw-col-ppb vs grid-nearest-ppb:")
for _,r in md.sample(5,random_state=1).iterrows():
    c=int(r['MUSICA0_colIndex']); lon360=r['lon']+360 if r['lon']<0 else r['lon']
    g=float(gi.sel(lat=r['lat'],lon=lon360,method='nearest'))
    print(f"  {r['AQS_code']:>16} raw={raw[18,c]:6.1f} grid={g:6.1f} d={g-raw[18,c]:+.1f}")

# ---- plots (daily mean) ----
try:
    import cartopy.crs as ccrs, cartopy.feature as cfeat
    HAVE_CARTO=True
except Exception as e:
    HAVE_CARTO=False; print("no cartopy:",e)
lon_p=lon2d-360.0  # all >180 -> -150..-50
ext=[-125,-66,23,50]
def basemap(ax):
    if HAVE_CARTO:
        ax.add_feature(cfeat.COASTLINE,lw=.5); ax.add_feature(cfeat.STATES,lw=.3,edgecolor='gray')
        ax.add_feature(cfeat.BORDERS,lw=.4); ax.set_extent(ext,ccrs.PlateCarree())
    else:
        ax.set_xlim(ext[0],ext[1]); ax.set_ylim(ext[2],ext[3])
proj=dict(projection=ccrs.PlateCarree()) if HAVE_CARTO else {}

base_dm=reg['BASE'].mean(0)
# 1) BASE map
fig,ax=plt.subplots(figsize=(8,5),subplot_kw=proj)
pc=ax.pcolormesh(lon_p,lat2d,base_dm,cmap='YlOrRd',vmin=20,vmax=60,
                 shading='auto',**(dict(transform=ccrs.PlateCarree()) if HAVE_CARTO else {}))
basemap(ax); plt.colorbar(pc,ax=ax,shrink=.8,label='surface O3 (ppb)')
ax.set_title(f'BASE daily-mean surface O3  {DAY}  (ne30 -> 1deg)')
plt.tight_layout(); plt.savefig(FIG+f'BASE_dailymean_{DAY}.png',dpi=130); plt.close()

# 2&3) differences
for k in ['noAnthro','noBB']:
    diff=base_dm-reg[k].mean(0)
    vmax=np.nanpercentile(np.abs(diff),99)
    fig,ax=plt.subplots(figsize=(8,5),subplot_kw=proj)
    pc=ax.pcolormesh(lon_p,lat2d,diff,cmap='RdBu_r',vmin=-vmax,vmax=vmax,
                     shading='auto',**(dict(transform=ccrs.PlateCarree()) if HAVE_CARTO else {}))
    basemap(ax); plt.colorbar(pc,ax=ax,shrink=.8,label='ΔO3 (ppb)')
    ax.set_title(f'BASE - {k}  daily-mean surface O3  {DAY}')
    plt.tight_layout(); plt.savefig(FIG+f'BASE_minus_{k}_{DAY}.png',dpi=130); plt.close()
    print(f"BASE-{k}: max +{np.nanmax(diff):.1f} / min {np.nanmin(diff):.1f} ppb")

print("\nfigures ->",FIG)
