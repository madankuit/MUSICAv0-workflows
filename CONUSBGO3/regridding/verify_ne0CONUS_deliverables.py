#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify_ne0CONUS_deliverables.py

Checks the ne0CONUSne30x8 CONUSBGO3 deliverables (BASE, noBBCONUS; Apr-Oct 2023 +
2024) before they are shared, and draws a comparison figure.

1. Merged hourly series: BASE and noBBCONUS have the same time axis each year,
   starting Apr 1 01:00 (end-of-hour stamps), and come from the intended cases.
2. Unified MDA8 files (1x1 and 0.15): scenarios, 428 dates, case names.
3. Point files: case/grid attributes; gridded vs point MDA8 (nearest cell) for
   both targets.
4. 0.15 vs 1x1: 0.15 seasonal means block-averaged to 1x1 against the 1x1 file.
5. Fire signal: BASE - noBBCONUS seasonal mean per year; the ne30 2023
   BASE - noBBCONUS on the same 1x1 cells for reference.

Figure -> BGO3_FIGURE_DIR/ne0CONUS/ne0CONUS_BASE_minus_noBBCONUS_MDA8_<year>.png

    python verify_ne0CONUS_deliverables.py <cdate YYYYMMDD>

    9 Oct 2026: VERSION 1.0
"""
import sys, pathlib, glob
import numpy as np, pandas as pd, xarray as xr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import cartopy.crs as ccrs, cartopy.feature as cfeat

_ROOT = next(p for p in pathlib.Path(__file__).resolve().parents
             if (p / 'config' / 'paths.py').exists())
sys.path.insert(0, str(_ROOT))
import config  # noqa: F401
from config.paths import (BGO3_EXPERIMENTS, BGO3_FIGURE_DIR, bgo3_merged_surfo3_path,
                          bgo3_unified_mda8_latest, ensure_dir)

CDATE = sys.argv[1]
E = BGO3_EXPERIMENTS['ne0CONUS']
FIG = ensure_dir(BGO3_FIGURE_DIR / 'ne0CONUS')

# 1. merged series
for yr in E['years']:
    t = {}
    for sc in E['scenarios']:
        p = bgo3_merged_surfo3_path(E['cases'][(sc, yr)], yr, E['merged_dir'])
        ds = xr.open_dataset(p)
        assert ds.O3.attrs['case_name'] == E['cases'][(sc, yr)] and ds.sizes['ncol'] == E['ncol']
        t[sc] = pd.DatetimeIndex(ds.time.values)
        print(f"1. {sc}{yr}: {pathlib.Path(p).name}  {t[sc][0]} .. {t[sc][-1]}  n={len(t[sc])}")
    common = t['BASE'].intersection(t['noBBCONUS'])
    assert t['BASE'][0] == t['noBBCONUS'][0] == pd.Timestamp(f'{yr}-04-01 01:00'), "different first stamp"
    assert len(common) >= 24 * 214, "time axes do not overlap over the season"
    print(f"   {yr}: common hourly stamps {len(common)}; BASE-only {len(t['BASE'].difference(common))}, "
          f"noBBCONUS-only {len(t['noBBCONUS'].difference(common))} (end-of-run tail only)")

# 2. unified files
U = {}
for tgt, (_, _, tag) in E['targets'].items():
    f = E['regrid_dir'] / (E['unified_stem'].format(tag=tag) + f'_c{CDATE}.nc')
    U[tgt] = xr.open_dataset(f)
    m = U[tgt]
    assert list(m.scenario.values) == E['scenarios'] and m.sizes['time'] == 428
    for (sc, yr), cn in E['cases'].items():
        assert f"{E['label_prefix']}{sc}{yr}: {cn}" in m.attrs['case_names']
    print(f"2. {f.name}: dims {dict(m.MDA8O3.sizes)}  NaN% {100*float(np.isnan(m.MDA8O3).mean()):.2f}")

# 3. points + grid-vs-point
P = {}
for (sc, yr), cn in E['cases'].items():
    f = glob.glob(str(E['points_dir'] / f"LocalTimeMDA8O3.{E['label_prefix']}{sc}{yr}.GivenMonitors.*.nc"))
    assert len(f) == 1, f
    P[(sc, yr)] = xr.open_dataset(f[0])
    assert P[(sc, yr)].attrs['case_name'] == cn and P[(sc, yr)].attrs['grid'] == E['grid']
print(f"3. point files: {len(P)}; sites per file {P[('BASE', 2023)].sizes['AQS_code']}")
for tgt, m in U.items():
    for (sc, yr), p in P.items():
        inbox = (p.lat >= 24) & (p.lat <= 50) & (p.lon >= -125) & (p.lon <= -66)
        p = p.isel(AQS_code=np.where(inbox.values)[0])
        g = m.MDA8O3.sel(scenario=sc, time=str(yr)).sel(
            lat=xr.DataArray(p.lat.values, dims='AQS_code'),
            lon=xr.DataArray(p.lon.values, dims='AQS_code'), method='nearest')
        pv = p.MDA8O3.transpose('AQS_code', 'time').values; gv = g.transpose('AQS_code', 'time').values
        ok = np.isfinite(pv) & np.isfinite(gv)
        print(f"   {tgt:4s} {sc}{yr}: gridded vs point r={np.corrcoef(pv[ok], gv[ok])[0, 1]:.3f} "
              f"bias={np.mean(gv[ok] - pv[ok]):+.2f} ppb (n sites {p.sizes['AQS_code']})")

# 4. 0.15 block-averaged to 1x1 vs 1x1 (seasonal means; MDA8 is nonlinear so not exact)
a1 = U['1x1'].MDA8O3.mean('time'); a15 = U['0.15'].MDA8O3.mean('time')
LA, LO = np.meshgrid(a15.lat.values, a15.lon.values, indexing='ij')
iy = np.round(LA).astype(int); ix = np.round(LO).astype(int)          # 1x1 cell whose centre is nearest
for sc in E['scenarios']:
    blk = pd.DataFrame({'lat': iy.ravel(), 'lon': ix.ravel(), 'v': a15.sel(scenario=sc).values.ravel()}) \
            .groupby(['lat', 'lon']).v.mean()
    x = a1.sel(scenario=sc).to_series().rename_axis(['lat', 'lon'])
    x.index = pd.MultiIndex.from_arrays([x.index.get_level_values(0).astype(int), x.index.get_level_values(1).astype(int)])
    j2 = pd.concat([x.rename('one'), blk.rename('blk')], axis=1, join='inner').dropna()
    print(f"4. {sc}: 0.15 block-mean vs 1x1 seasonal mean  r={np.corrcoef(j2.one, j2.blk)[0, 1]:.3f}  "
          f"bias={np.mean(j2.blk - j2.one):+.2f} ppb  (n={len(j2)} cells)")

# 5. fire signal
ne30 = xr.open_dataset(bgo3_unified_mda8_latest()).MDA8O3
for yr in E['years']:
    for tgt, m in U.items():
        d = (m.MDA8O3.sel(scenario='BASE', time=str(yr)) - m.MDA8O3.sel(scenario='noBBCONUS', time=str(yr))).mean('time')
        print(f"5. {tgt:4s} {yr}: BASE-noBBCONUS seasonal mean: CONUS mean {float(d.mean()):.2f}, "
              f"max {float(d.max()):.2f}, min {float(d.min()):.2f} ppb")
d30 = (ne30.sel(scenario='BASE', time='2023') - ne30.sel(scenario='noBBCONUS', time='2023')).mean('time')
d0 = (U['1x1'].MDA8O3.sel(scenario='BASE', time='2023') - U['1x1'].MDA8O3.sel(scenario='noBBCONUS', time='2023')).mean('time')
ok = np.isfinite(d30.values) & np.isfinite(d0.values)
print(f"   2023 1x1 BASE-noBBCONUS: ne30 CONUS mean {float(d30.mean()):.2f} vs ne0CONUS {float(d0.mean()):.2f} ppb; "
      f"pattern r={np.corrcoef(d30.values[ok], d0.values[ok])[0, 1]:.2f} (different model setups)")

# figure: per year, BASE and BASE-noBBCONUS on 1x1 and 0.15 (+ ne30 1x1 for 2023)
for yr in E['years']:
    panels = [(f'ne0CONUS 1x1', U['1x1']), ('ne0CONUS 0.15', U['0.15'])]
    fig, axs = plt.subplots(1, 3 if yr == 2023 else 2, figsize=(15 if yr == 2023 else 10.5, 3.8),
                            subplot_kw={'projection': ccrs.PlateCarree()}, constrained_layout=True)
    fields = [(n, (m.MDA8O3.sel(scenario='BASE', time=str(yr)) - m.MDA8O3.sel(scenario='noBBCONUS', time=str(yr))).mean('time'))
              for n, m in panels]
    if yr == 2023:
        fields.append(('ne30 1x1 (CAMS/MERRA-2 setup)', d30))
    v = float(np.ceil(max(np.nanpercentile(f.values, 99.5) for _, f in fields)))
    for ax, (n, f) in zip(np.atleast_1d(axs), fields):
        pc = ax.pcolormesh(f.lon, f.lat, f, cmap='RdBu_r', vmin=-v, vmax=v, shading='auto',
                           transform=ccrs.PlateCarree())
        ax.add_feature(cfeat.COASTLINE, lw=.5); ax.add_feature(cfeat.STATES, lw=.3, edgecolor='gray')
        ax.add_feature(cfeat.BORDERS, lw=.4); ax.set_extent([-125, -66, 24, 50], ccrs.PlateCarree())
        ax.set_title(n, fontsize=11, fontweight='bold')
    fig.colorbar(pc, ax=axs, orientation='horizontal', shrink=.5, pad=.03, extend='both') \
       .set_label('BASE $-$ noBBCONUS seasonal-mean MDA8 O$_3$ (ppb)')
    fig.suptitle(f'US-fire MDA8 O$_3$ (BASE $-$ noBBCONUS), Apr$-$Oct {yr}', fontsize=13, fontweight='bold')
    out = FIG / f'ne0CONUS_BASE_minus_noBBCONUS_MDA8_{yr}.png'
    fig.savefig(out, dpi=150); plt.close(fig); print('saved', out)
print('all checks passed')
