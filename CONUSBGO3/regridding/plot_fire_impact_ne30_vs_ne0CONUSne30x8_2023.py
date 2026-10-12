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
    ('ne30np4 (~111 km global)\nCAMS-GLOB-ANT v6.2 + MERRA-2 met\n1° output', bgo3_unified_mda8_latest()),
    ('ne0CONUSne30x8 (~14 km over CONUS)\nNEI2022v2 hourly (CONUS) + CAMS v6.2 + GEOS-5 met\n0.15° output',
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

fig, axs = plt.subplots(2, 2, figsize=(13.5, 7.4), subplot_kw={'projection': ccrs.PlateCarree()},
                        constrained_layout=True)
for r, (lab, f, base, fire) in enumerate(fields):
    for c, (da, cmap, lo, hi, ttl) in enumerate([
            (base, 'YlOrRd', bmin, bmax, 'BASE MDA8 O$_3$'),
            (fire, 'RdBu_r', -vf, vf, 'Fire impact: BASE $-$ noBBCONUS')]):
        ax = axs[r, c]
        pc = ax.pcolormesh(da.lon, da.lat, da, cmap=cmap, vmin=lo, vmax=hi, shading='auto',
                           transform=ccrs.PlateCarree())
        ax.add_feature(cfeat.COASTLINE, lw=.5); ax.add_feature(cfeat.STATES, lw=.3, edgecolor='gray')
        ax.add_feature(cfeat.BORDERS, lw=.4); ax.set_extent(EXT, ccrs.PlateCarree())
        if r == 0:
            ax.set_title(ttl, fontsize=12, fontweight='bold')
        if c == 0:
            ax.text(-0.03, 0.5, lab, transform=ax.transAxes, rotation=90, ha='right', va='center',
                    fontsize=9.5, fontweight='bold')
        ax.text(0.01, 0.03, f'box mean {boxmean(da):.2f} ppb\nmax {float(da.max()):.1f} ppb',
                transform=ax.transAxes, fontsize=8.5, va='bottom',
                bbox=dict(fc='white', ec='0.6', alpha=.85))
        if r == 1:
            fig.colorbar(pc, ax=axs[:, c], orientation='horizontal', shrink=.8, pad=.02, aspect=40,
                         extend='both').set_label('seasonal-mean MDA8 O$_3$ (ppb)' if c == 0 else
                                                  '$\\Delta$MDA8 O$_3$ (ppb)')
fig.suptitle(f'US-fire impact on surface O$_3$, Apr$-$Oct {YEAR}: ne30np4 vs ne0CONUSne30x8',
             fontsize=13, fontweight='bold')
out = FIG / f'FireImpact_BASEminusNoBBCONUS_ne30np4_vs_ne0CONUSne30x8_AprOct{YEAR}.png'
fig.savefig(out, dpi=170); plt.close(fig); print('saved', out)
