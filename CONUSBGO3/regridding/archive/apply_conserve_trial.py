#!/usr/bin/env python3
"""
apply_conserve_trial.py  (ARCHIVE - method development, 8 Jul 2026)

Single-day trial of the ne30np4 -> 1x1 **mass-conservative** regrid, run before
committing to the full Apr-Oct 2022-2023 production job. It:

  1. builds the sparse conservative operator from the offline ESMF weights,
  2. regrids one day (2022-07-20) of merged surface O3 for BASE/noAnthro/noBB,
  3. cross-checks regridded cell values against the raw ne30 column values at
     a sample of monitors - the check that established the method was sound,
  4. plots BASE plus the two scenario differences for that day.

This is the approach that was adopted; the production version is
../Regrid_ne30_surfO3_to_1x1_conserve.py. Kept as the record of the validation
that justified it, and as a quick smoke test if the weight file is ever
regenerated. Compare with trial_regrid_linear.py, the rejected alternative.

MODIFICATION HISTORY:
    8 Jul 2026: VERSION 1.0
    31 Aug 2026: VERSION 1.1 - paths moved to config/paths.py
"""
import numpy as np, pandas as pd, xarray as xr
from scipy.sparse import coo_matrix
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

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
    WEIGHTS_NE30_TO_1X1,
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

# --- destination 1x1 axes (must match FV1x1grid_info) ---
lat2d=np.arange(-90,90+1,1.0); lon2d=np.arange(0,360+1,1.0)   # 181, 361
nlat,nlon=lat2d.size,lon2d.size; ndst=nlat*nlon

# --- build sparse conservative operator M (ndst x nsrc) ---
w=xr.open_dataset(str(WEIGHTS_NE30_TO_1X1)); col=w["col"].values-1; row=w["row"].values-1; S=w["S"].values
nsrc=48602
M=coo_matrix((S,(row,col)),shape=(ndst,nsrc)).tocsr()
wsum=np.asarray(M.sum(axis=1)).ravel()          # ~1 where mapped, 0 where unmapped
valid=wsum>0.5

# --- locate each dest index on the lat/lon axes via area-weighted src coords ---
slat=xr.open_dataset(h2ex)['lat'].values; slon=xr.open_dataset(h2ex)['lon'].values
dlat=M.dot(slat); dlon=M.dot(slon)
ilat=np.full(ndst,-1); ilon=np.full(ndst,-1)
ilat[valid]=np.rint(dlat[valid]+90).astype(int)
ilon[valid]=np.rint(dlon[valid]).astype(int)
ok=valid&(ilat>=0)&(ilat<nlat)&(ilon>=0)&(ilon<nlon)

def to_grid(vec_src):
    d=M.dot(vec_src); g=np.full((nlat,nlon),np.nan)
    g[ilat[ok],ilon[ok]]=d[ok]/wsum[ok]
    return g

# --- regrid daily-mean for each case ---
reg={}
for k,p in merged.items():
    o=xr.open_dataset(p).sel(time=slice(DAY+"T00",DAY+"T23"))['O3'].values*1e9  # (24,ncol) ppb
    reg[k]=to_grid(o.mean(0))
    print(f"{k:9s} dailymean ppb: min={np.nanmin(reg[k]):.1f} max={np.nanmax(reg[k]):.1f}")

# --- cross-check BASE grid cell vs raw column at monitors ---
md=pd.read_csv(BGO3_MONITOR_COLIDX)
md=md[md['MUSICA0_colIndex']!='Find None'].dropna(subset=['MUSICA0_colIndex'])
md['MUSICA0_colIndex']=md['MUSICA0_colIndex'].astype(int)
raw=xr.open_dataset(merged['BASE']).sel(time=slice(DAY+"T00",DAY+"T23"))['O3'].values*1e9
base_h18=to_grid(raw[18])
gi=xr.DataArray(base_h18,dims=['lat','lon'],coords={'lat':lat2d,'lon':lon2d})
print("\ncross-check hour18  raw-col-ppb vs 1x1-CONSERVE-cell-ppb:")
for _,r in md.sample(5,random_state=1).iterrows():
    c=int(r['MUSICA0_colIndex']); lon360=r['lon']+360 if r['lon']<0 else r['lon']
    g=float(gi.sel(lat=r['lat'],lon=lon360,method='nearest'))
    print(f"  {r['AQS_code']:>16} raw={raw[18,c]:6.1f} cell={g:6.1f} d={g-raw[18,c]:+.1f}")

# --- plots ---
import cartopy.crs as ccrs, cartopy.feature as cfeat
lon_p=np.where(lon2d>180,lon2d-360,lon2d)
order=np.argsort(lon_p); lon_ps=lon_p[order]
ext=[-125,-66,23,50]
def basemap(ax):
    ax.add_feature(cfeat.COASTLINE,lw=.5); ax.add_feature(cfeat.STATES,lw=.3,edgecolor='gray')
    ax.add_feature(cfeat.BORDERS,lw=.4); ax.set_extent(ext,ccrs.PlateCarree())
proj=dict(projection=ccrs.PlateCarree())
base_dm=reg['BASE']
fig,ax=plt.subplots(figsize=(8,5),subplot_kw=proj)
pc=ax.pcolormesh(lon_ps,lat2d,base_dm[:,order],cmap='YlOrRd',vmin=20,vmax=60,shading='auto',transform=ccrs.PlateCarree())
basemap(ax); plt.colorbar(pc,ax=ax,shrink=.8,label='surface O3 (ppb)')
ax.set_title(f'BASE daily-mean surface O3  {DAY}  (ne30 -> 1x1, CONSERVATIVE)')
plt.tight_layout(); plt.savefig(FIG+f'CONSERVE_BASE_dailymean_{DAY}.png',dpi=130); plt.close()
for k in ['noAnthro','noBB']:
    diff=base_dm-reg[k]; vmax=np.nanpercentile(np.abs(diff),99)
    fig,ax=plt.subplots(figsize=(8,5),subplot_kw=proj)
    pc=ax.pcolormesh(lon_ps,lat2d,diff[:,order],cmap='RdBu_r',vmin=-vmax,vmax=vmax,shading='auto',transform=ccrs.PlateCarree())
    basemap(ax); plt.colorbar(pc,ax=ax,shrink=.8,label='ΔO3 (ppb)')
    ax.set_title(f'BASE - {k}  daily-mean surface O3  {DAY}  (CONSERVATIVE)')
    plt.tight_layout(); plt.savefig(FIG+f'CONSERVE_BASE_minus_{k}_{DAY}.png',dpi=130); plt.close()
    print(f"BASE-{k}: max +{np.nanmax(diff):.1f} / min {np.nanmin(diff):.1f} ppb")
print("\nfigures ->",FIG)
