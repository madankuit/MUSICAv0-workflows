#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compare_unified_versions.py

Compares two versions of each unified MDA8 file (e.g. c20261010 vs c20261012, the
per-year time-zone fix). For every (experiment, target, scenario, year): cells whose
MDA8 changed by more than 1e-3 ppb, and whether EVERY changed cell is one whose
summertime UTC offset differs between that year's rules and 2022's (i.e. Mexico).
Also confirms no monitor's offset depends on the year (point files unaffected).

    python compare_unified_versions.py <old cdate> <new cdate>

    12 Oct 2026: VERSION 1.0
"""
import sys, pathlib, glob
import numpy as np, pandas as pd, xarray as xr
from timezonefinder import TimezoneFinder
import pytz
from datetime import datetime
_ROOT = next(p for p in pathlib.Path(__file__).resolve().parents if (p / 'config' / 'paths.py').exists())
sys.path.insert(0, str(_ROOT))
import config  # noqa: F401
from config.paths import BGO3_EXPERIMENTS

OLD, NEW = sys.argv[1], sys.argv[2]
tf = TimezoneFinder()
def off(tz, lon, yr):
    if tz is None:
        return int(round(lon / 15.0))
    return int(pytz.timezone(tz).localize(datetime(yr, 7, 1), is_dst=True).utcoffset().total_seconds() // 3600)
ok = True
for ex, E in BGO3_EXPERIMENTS.items():
    for tgt, (_, _, tag) in E['targets'].items():
        stem = E['unified_stem'].format(tag=tag)
        fo = glob.glob(str(E['regrid_dir'] / f'{stem}_c{OLD}.nc')); fn = glob.glob(str(E['regrid_dir'] / f'{stem}_c{NEW}.nc'))
        if not (fo and fn):
            print(f'{ex} {tgt}: missing {OLD if not fo else NEW}'); ok = False; continue
        o = xr.open_dataset(fo[0]); n = xr.open_dataset(fn[0])
        o = o.assign_coords(scenario=[{'noBB': 'noBBCONUS'}.get(str(v), str(v)) for v in o.scenario.values])
        assert list(o.scenario.values) == list(n.scenario.values) and np.array_equal(o.time, n.time)
        LA, LO = np.meshgrid(n.lat.values, n.lon.values, indexing='ij')
        tzs = np.array([tf.timezone_at(lat=a, lng=b) for a, b in zip(LA.ravel(), LO.ravel())], dtype=object).reshape(LA.shape)
        for yr in E['years']:
            moved = np.array([[off(tzs[i, j], LO[i, j], yr) != off(tzs[i, j], LO[i, j], 2022)
                               for j in range(LA.shape[1])] for i in range(LA.shape[0])])
            for sc in n.scenario.values:
                a = o.MDA8O3.sel(scenario=sc, time=str(yr)).values; b = n.MDA8O3.sel(scenario=sc, time=str(yr)).values
                ch = (np.abs(a - b) > 1e-3).any(0) | (np.isnan(a) != np.isnan(b)).any(0)
                bad = ch & ~moved
                zones = sorted(set(tzs[ch & moved])) if (ch & moved).any() else []
                print(f'{ex} {tgt} {sc}{yr}: changed cells {int(ch.sum())} (all tz-offset-moved: {not bad.any()}); '
                      f'offset-moved cells {int(moved.sum())}; zones {zones[:6]}; max |d| {np.nanmax(np.abs(a - b)):.2f} ppb')
                ok &= not bad.any()
for ex, E in BGO3_EXPERIMENTS.items():
    m = pd.read_csv(E.get('colidx_exact', E['colidx']))
    var = [(r.AQS_code, sorted({off(tf.timezone_at(lat=r.lat, lng=r.lon), r.lon, y) for y in range(2022, 2027)}))
           for r in m.itertuples()]
    multi = [v for v in var if len(v[1]) > 1]
    print(f'{ex} monitors with a year-dependent offset (2022-2026): {len(multi)} {multi[:5]}')
    ok &= not multi
print('VERSIONS CONSISTENT' if ok else 'UNEXPECTED DIFFERENCES')
sys.exit(0 if ok else 1)
