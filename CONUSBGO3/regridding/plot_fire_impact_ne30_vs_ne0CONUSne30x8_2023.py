#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_fire_impact_ne30_vs_ne0CONUSne30x8_2023.py

US-fire impact on surface O3, Apr-Oct 2023, in the two CONUSBGO3 model setups:

  row 1  ne30np4 (~111 km)             CAMS-GLOB-ANT v6.2 + MERRA-2 met   (1x1 deliverable)
  row 2  ne0CONUSne30x8 (~14 km CONUS) NEI2022v2 + CAMS v6.2 + GEOS-5 met (0.15x0.15 deliverable)

  col 1  BASE seasonal-mean MDA8 O3              (shared colour scale)
  col 2  fire impact = BASE - noBBCONUS          (shared diverging scale)

Each panel is annotated with its cos(lat)-weighted mean over the CONUS box
(24-50N, 125-66W, land and ocean). Reads the newest unified MDA8 file of each
experiment; nothing is hard-coded except the year.

    python plot_fire_impact_ne30_vs_ne0CONUSne30x8_2023.py

    12 Oct 2026: VERSION 1.0
"""
import sys, pathlib, glob, re
import numpy as np, xarray as xr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import cartopy.crs as ccrs, cartopy.feature as cfeat

_ROOT = next(p for p in pathlib.Path(__file__).resolve().parents
             if (p / 'config' / 'paths.py').exists())
sys.path.insert(0, str(_ROOT))
import config  # noqa: F401
from config.paths import BGO3_EXPERIMENTS, BGO3_FIGURE_DIR, bgo3_unified_mda8_latest, ensure_dir

YEAR = 2023
EXT = [-125, -66, 24, 50]
FIG = ensure_dir(BGO3_FIGURE_DIR / 'ne0CONUSne30x8')


def newest(exp, target):
    E = BGO3_EXPERIMENTS[exp]
    stem = E['unified_stem'].format(tag=E['targets'][target][2])
    hits = glob.glob(str(E['regrid_dir'] / f'{stem}_c*.nc'))
    return max(hits, key=lambda h: re.search(r'_c(\d{8})\.nc$', h).group(1))


ROWS = [
    ('ne30np4 (~111 km global)  |  CAMS-GLOB-ANT v6.2 anthro + MERRA-2 met  |  1° output',
     bgo3_unified_mda8_latest()),
    ('ne0CONUSne30x8 (~14 km over CONUS)  |  NEI2022v2 hourly (CONUS) + CAMS v6.2 anthro + GEOS-5 met  |  0.15° output',
     newest('ne0CONUSne30x8', '0.15')),
]

fields = []
for lab, f in ROWS:
    m = xr.open_dataset(f).MDA8O3.sel(time=str(YEAR))
    assert m.sizes['time'] == 214, (f, m.sizes['time'])
    base = m.sel(scenario='BASE').mean('time')
    fire = base - m.sel(scenario='noBBCONUS').mean('time')
    fields.append((lab, f, base, fire))
    print(f"{pathlib.Path(f).name}: case_names -> "
          f"{[c.strip() for c in xr.open_dataset(f).attrs['case_names'].split(';') if str(YEAR) in c]}")


def boxmean(da):
    w = np.cos(np.deg2rad(da.lat))
    return float(da.weighted(w).mean())


allb = np.concatenate([b.values.ravel() for _, _, b, _ in fields])
bmin, bmax = np.floor(np.nanpercentile(allb, 2)), np.ceil(np.nanpercentile(allb, 98))
vf = float(np.ceil(max(np.nanpercentile(np.abs(d.values), 99.5) for *_, d in fields)))

fig = plt.figure(figsize=(14, 8.4))
X0, W = [0.03, 0.515], 0.455                    # panel columns
Y0, H = [0.555, 0.15], 0.335                    # panel rows (top, bottom)
cbars = {}
for r, (lab, f, base, fire) in enumerate(fields):
    fig.text(0.5, Y0[r] + H + 0.022, lab, ha='center', va='bottom', fontsize=11.5, fontweight='bold',
             color='#1f3b73' if r == 0 else '#7a1f1f')
    for c, (da, cmap, lo, hi, ttl) in enumerate([
            (base, 'YlOrRd', bmin, bmax, 'BASE MDA8 O$_3$'),
            (fire, 'RdBu_r', -vf, vf, 'Fire impact: BASE $-$ noBBCONUS')]):
        ax = fig.add_axes([X0[c], Y0[r], W, H], projection=ccrs.PlateCarree())
        pc = ax.pcolormesh(da.lon, da.lat, da, cmap=cmap, vmin=lo, vmax=hi, shading='auto',
                           transform=ccrs.PlateCarree())
        ax.add_feature(cfeat.COASTLINE, lw=.5); ax.add_feature(cfeat.STATES, lw=.3, edgecolor='gray')
        ax.add_feature(cfeat.BORDERS, lw=.4); ax.set_extent(EXT, ccrs.PlateCarree())
        ax.set_title(ttl, fontsize=10.5)
        ax.text(0.01, 0.03, f'box mean {boxmean(da):.2f} ppb\nmax {float(da.max()):.1f} ppb',
                transform=ax.transAxes, fontsize=8.5, va='bottom',
                bbox=dict(fc='white', ec='0.6', alpha=.85))
        cbars[c] = pc
for c, lab in [(0, 'seasonal-mean MDA8 O$_3$ (ppb)'), (1, '$\\Delta$MDA8 O$_3$ (ppb)')]:
    cax = fig.add_axes([X0[c] + 0.06, 0.075, W - 0.12, 0.022])
    fig.colorbar(cbars[c], cax=cax, orientation='horizontal', extend='both').set_label(lab)
fig.suptitle(f'US-fire impact on surface O$_3$, Apr$-$Oct {YEAR}: ne30np4 vs ne0CONUSne30x8',
             fontsize=14, fontweight='bold', y=0.985)
out = FIG / f'FireImpact_BASEminusNoBBCONUS_ne30np4_vs_ne0CONUSne30x8_AprOct{YEAR}.png'
fig.savefig(out, dpi=170); plt.close(fig); print('saved', out)
