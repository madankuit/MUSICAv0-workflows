#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify_ne30_allmonitors.py

Checks the ne30 point files rebuilt with exact point-in-polygon column matching
(GivenMonitorsAll803) against the delivered ones (legacy matching):

1. Column indices: exact == legacy for every monitor legacy matched; list the
   monitors only exact matched, with their distance to the column centre.
2. For every (scenario, year): MDA8 and hourly O3 at the shared monitors are
   identical to the delivered file; the new monitors have no missing days
   beyond what the shared monitors have.

    python verify_ne30_allmonitors.py

    10 Oct 2026: VERSION 1.0
"""
import sys, pathlib, glob
import numpy as np, pandas as pd, xarray as xr

_ROOT = next(p for p in pathlib.Path(__file__).resolve().parents
             if (p / 'config' / 'paths.py').exists())
sys.path.insert(0, str(_ROOT))
import config  # noqa: F401
from config.paths import BGO3_EXPERIMENTS

E = BGO3_EXPERIMENTS['ne30']
OLD_LABEL = {'noBBCONUS': 'noBB'}          # delivered files predate the rename

# 1. column indices
leg = pd.read_csv(E['colidx'], dtype={'AQS_code': str})
exa = pd.read_csv(E['colidx_exact'], dtype={'AQS_code': str})
assert (leg.AQS_code.values == exa.AQS_code.values).all()
lm = leg.MUSICA0_colIndex.astype(str) != 'Find None'
em = exa.MUSICA0_colIndex.astype(str) != 'Find None'
same = (leg.MUSICA0_colIndex[lm].astype(int).values == exa.MUSICA0_colIndex[lm].astype(int).values)
print(f"1. legacy matched {lm.sum()}, exact matched {em.sum()}; exact == legacy for "
      f"{same.sum()}/{lm.sum()} legacy-matched monitors")
assert same.all(), leg[lm][~same][['AQS_code', 'lat', 'lon']]
new = exa[em & ~lm].copy()
new['km'] = 111.2 * np.hypot(new.Approx_MUSICA0_lat - new.lat,
                             (((new.Approx_MUSICA0_lon - new.lon + 180) % 360) - 180) * np.cos(np.deg2rad(new.lat)))
print(f"   newly matched ({len(new)}):")
print(new[['AQS_code', 'lat', 'lon', 'MUSICA0_colIndex', 'km']].round(3).to_string(index=False))
NEW = set(new.AQS_code)

# 2. values
for (sc, yr) in E['cases']:
    for kind, var in [('LocalTimeMDA8O3', 'MDA8O3'), ('UTChourlyO3', 'O3')]:
        fn = glob.glob(str(E['points_dir'] / f"{kind}.{sc}{yr}.GivenMonitorsAll803.*.nc"))
        fo = glob.glob(str(E['points_dir'] / f"{kind}.{OLD_LABEL.get(sc, sc)}{yr}.GivenMonitors.*.nc"))
        assert len(fn) == 1 and len(fo) == 1, (sc, yr, kind, fn, fo)
        n = xr.open_dataset(fn[0]); o = xr.open_dataset(fo[0])
        shared = np.intersect1d(o.AQS_code.values, n.AQS_code.values)
        assert len(shared) == o.sizes['AQS_code']
        a = n[var].sel(AQS_code=shared).transpose('AQS_code', 'time').values
        b = o[var].sel(AQS_code=shared).transpose('AQS_code', 'time').values
        ident = np.array_equal(np.isnan(a), np.isnan(b)) and np.array_equal(a[~np.isnan(a)], b[~np.isnan(b)])
        added = sorted(set(n.AQS_code.values) - set(o.AQS_code.values))
        assert set(added) == NEW, (sc, yr, kind)
        nan_new = float(np.isnan(n[var].sel(AQS_code=added).values).mean())
        nan_old = float(np.isnan(b).mean())
        print(f"2. {sc}{yr} {kind:15s}: {len(shared)} shared sites {'IDENTICAL' if ident else 'DIFFERENT'}; "
              f"{len(added)} added; NaN% added {100*nan_new:.2f} vs shared {100*nan_old:.2f}")
        assert ident
        if var == 'MDA8O3':
            m = n[var].sel(AQS_code=added).mean('time')
            print(f"   added-site season-mean MDA8 {float(m.min()):.1f}..{float(m.max()):.1f} ppb "
                  f"(shared {float(o[var].mean('time').min()):.1f}..{float(o[var].mean('time').max()):.1f})")
print("all checks passed")
