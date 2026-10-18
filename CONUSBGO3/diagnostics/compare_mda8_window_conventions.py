#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compare_mda8_window_conventions.py

How much does the MDA8 window choice matter? The pipeline (since Sep 2025) takes
the daily maximum over the 17 8-h windows whose LAST hourly value is labelled
07..23 local, i.e. windows covering 23:00(prev)-07:00 ... 15:00-23:00. The EPA
2015 definition (40 CFR 50 App. U) uses the 17 windows BEGINNING 07:00..23:00
(07:00-15:00 ... 23:00-07:00 next day); with end-of-hour stamps those are the
windows whose last value is labelled 15..31 h after local midnight.

From the 1x1 hourly files: recompute MDA8 both ways (>= 6 valid hours per window,
>= 13 of 17 windows per day), check the 'current' one reproduces the unified file,
and report EPA - current per scenario-year, plus the effect on seasonal-mean
BASE - noBBGlobal and BASE - noBBCONUS. Apr 1 - Oct 30 (the hourly files end
Oct 31 23 UTC, too early for the EPA windows of Oct 31). Writes a map figure.

    python compare_mda8_window_conventions.py [experiment]     # default ne30

    18 Oct 2026: VERSION 1.0
"""
import sys, pathlib, glob
import numpy as np, pandas as pd, xarray as xr
from timezonefinder import TimezoneFinder
import pytz
from datetime import datetime
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import cartopy.crs as ccrs, cartopy.feature as cfeat
_ROOT = next(p for p in pathlib.Path(__file__).resolve().parents if (p / 'config' / 'paths.py').exists())
sys.path.insert(0, str(_ROOT))
import config  # noqa: F401
from config.paths import BGO3_EXPERIMENTS, BGO3_FIGURE_DIR, ensure_dir

EX = sys.argv[1] if len(sys.argv) > 1 else 'ne30'
E = BGO3_EXPERIMENTS[EX]; tag = E['targets']['1x1'][2]; lp = E['label_prefix']
U = xr.open_dataset(sorted(glob.glob(str(E['regrid_dir'] / (E['unified_stem'].format(tag=tag) + '_c20261012.nc'))))[0]).MDA8O3
tf = TimezoneFinder()
LA, LO = np.meshgrid(U.lat.values, U.lon.values, indexing='ij'); lat, lon = LA.ravel(), LO.ravel()
tzn = [tf.timezone_at(lat=a, lng=b) for a, b in zip(lat, lon)]
def offsets(yr):
    return np.array([int(round(b / 15)) if z is None else
                     int(pytz.timezone(z).localize(datetime(yr, 7, 1), is_dst=True).utcoffset().total_seconds() // 3600)
                     for z, b in zip(tzn, lon)])

def mda8(H, t_utc, off, yr, end_hours):
    """H (nt, ncell) ppb -> (ndays, ncell). end_hours: labels of the windows' last value after local midnight."""
    days = pd.date_range(f'{yr}-04-01', f'{yr}-10-30')
    out = np.full((len(days), H.shape[1]), np.nan, dtype='float32')
    for o in np.unique(off):
        idx = np.where(off == o)[0]
        tl = pd.DatetimeIndex(t_utc) + pd.Timedelta(hours=int(o))
        grid = pd.date_range(days[0] - pd.Timedelta(hours=8), days[-1] + pd.Timedelta(hours=32), freq='h')
        pos = grid.get_indexer(tl); keep = pos >= 0
        A = np.full((len(grid), len(idx)), np.nan); A[pos[keep]] = H[keep][:, idx]
        v = np.nan_to_num(A); n = np.isfinite(A).astype(float)
        cs = np.vstack([np.zeros((1, len(idx))), np.cumsum(v, 0)]); cn = np.vstack([np.zeros((1, len(idx))), np.cumsum(n, 0)])
        S = cs[8:] - cs[:-8]; N = cn[8:] - cn[:-8]                 # window ending at grid[k+7]
        M = np.where(N >= 6, S / np.maximum(N, 1), np.nan)
        endt = grid[7:]
        for d, day in enumerate(days):
            k = endt.get_indexer(day + pd.to_timedelta(end_hours, 'h'))
            W = M[k]
            ok = np.isfinite(W).sum(0) >= 13
            with np.errstate(all='ignore'):
                out[d, idx[ok]] = np.nanmax(W[:, ok], 0)
    return days, out

CUR, EPA = list(range(7, 24)), list(range(15, 32))
res = {}
print(f'{EX}: EPA-2015 windows minus current windows, daily MDA8 at {lat.size} 1x1 cells, Apr 1 - Oct 30')
for (sc, yr), cn in E['cases'].items():
    h = sorted(glob.glob(str(E['hourly_dir'] / f'{tag}_UTChourlySurfO3.{lp}{sc}{yr}.*_c20261012.nc')))[0]
    ds = xr.open_dataset(h); t = ds.time.values; H = ds.O3.values.reshape(len(t), -1).astype('float64')
    off = offsets(yr)
    days, cur = mda8(H, t, off, yr, CUR); _, epa = mda8(H, t, off, yr, EPA)
    ref = U.sel(scenario=sc, time=slice(f'{yr}-04-01', f'{yr}-10-30')).values.reshape(len(days), -1)
    repro = float(np.nanmax(np.abs(cur - ref)))
    d = epa - cur; ad = np.abs(d)
    k = np.nanargmax(ad); di, ci = np.unravel_index(k, d.shape)
    res[(sc, yr)] = (cur, epa)
    print(f'  {lp}{sc}{yr}: reproduce-file max|d| {repro:.1e} | EPA-current mean {np.nanmean(d):+.2f} ppb, '
          f'|d|>1 ppb {100*np.nanmean(ad > 1):.1f}%, >5 ppb {100*np.nanmean(ad > 5):.2f}%, '
          f'unchanged {100*np.nanmean(ad < 1e-3):.1f}% of cell-days; max |d| {ad.flat[k]:.1f} ppb '
          f'({days[di].date()}, {lat[ci]:.0f}N {lon[ci]:.0f}E); seasonal-mean change '
          f'{np.nanmin(np.nanmean(d,0)):+.2f}..{np.nanmax(np.nanmean(d,0)):+.2f} ppb')

def sm(a): return np.nanmean(a, 0).reshape(LA.shape)
print('\nEffect on seasonal-mean scenario differences (1x1 cells, box mean / largest change in any cell):')
fig_rows = []
for other in [s for s in E['scenarios'] if s != 'BASE']:
    for yr in E['years']:
        c = sm(res[('BASE', yr)][0]) - sm(res[(other, yr)][0]); e = sm(res[('BASE', yr)][1]) - sm(res[(other, yr)][1])
        print(f'  BASE-{other} {yr}: current {np.nanmean(c):.2f}, EPA {np.nanmean(e):.2f} ppb (box mean); '
              f'max |change| {np.nanmax(np.abs(e - c)):.2f} ppb; pattern r {np.corrcoef(c.ravel(), e.ravel())[0,1]:.3f}')
        if other in ('noBBGlobal', 'noBBCONUS'):
            fig_rows.append((f'BASE$-${other} {yr}', c, e))

# figure: BASE seasonal mean change + BB-difference current vs EPA
if fig_rows:
    yr0 = E['years'][-1]
    db = sm(res[('BASE', yr0)][1]) - sm(res[('BASE', yr0)][0])
    r = [x for x in fig_rows if 'noBBGlobal' in x[0] and str(yr0) in x[0]] or fig_rows[-1:]
    lab, c, e = r[0]
    fig, axs = plt.subplots(1, 4, figsize=(20, 3.9), subplot_kw={'projection': ccrs.PlateCarree()}, constrained_layout=True)
    v1 = float(np.nanpercentile(np.abs(db), 99)); v2 = float(np.nanpercentile(np.abs(np.r_[c.ravel(), e.ravel()]), 99))
    for ax, f, ttl, cm, v in [(axs[0], db, f'BASE seasonal-mean MDA8 {yr0}: EPA $-$ current', 'RdBu_r', v1),
                              (axs[1], c, f'{lab}, current windows', 'RdBu_r', v2),
                              (axs[2], e, f'{lab}, EPA-2015 windows', 'RdBu_r', v2),
                              (axs[3], e - c, f'{lab}: EPA $-$ current', 'RdBu_r', max(float(np.nanpercentile(np.abs(e - c), 99)), 0.05))]:
        pc = ax.pcolormesh(U.lon, U.lat, f, cmap=cm, vmin=-v, vmax=v, shading='auto', transform=ccrs.PlateCarree())
        ax.add_feature(cfeat.COASTLINE, lw=.5); ax.add_feature(cfeat.STATES, lw=.3, edgecolor='gray')
        ax.set_extent([-125, -66, 24, 50], ccrs.PlateCarree()); ax.set_title(ttl, fontsize=10)
        fig.colorbar(pc, ax=ax, orientation='horizontal', shrink=.85, pad=.03).set_label('ppb')
    out = ensure_dir(BGO3_FIGURE_DIR / 'mda8_window_check') / f'MDA8_windows_EPA_vs_current_{EX}.png'
    fig.savefig(out, dpi=140); print('saved', out)
