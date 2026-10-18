#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
check_share_folder.py

Pre-send check of the folder shared with the GAM group (Oct 2026). Reads the files
themselves (h5netcdf backend, so it runs on a laptop without netCDF4):

1. Every global and variable attribute of every .nc, printed in full.
2. Unified MDA8: scenario list, 428 dates, 1x1 CONUS grid, units, no bare 'noBB',
   and wording that still says 'fire' instead of 'biomass burning'.
3. Drop-in compatibility with the July file: same variable, dimension order, dtype,
   coordinates and fill value; BASE / noAnthro identical; noBBCONUS == July noBB in
   2022; 2023 differences only in the listed cells.
4. Hourly files: layout, units, time span, case, scenario.

    python check_share_folder.py <share folder> <July unified MDA8 file>

    18 Oct 2026: VERSION 1.0
"""
import sys, re, glob, pathlib
import numpy as np, pandas as pd, xarray as xr

FOLDER = pathlib.Path(sys.argv[1]); JULY = pathlib.Path(sys.argv[2])
op = lambda p: xr.open_dataset(p, engine='h5netcdf')
BARE = re.compile(r'noBB(?![A-Za-z])')
probs = []
def chk(ok, msg):
    print(('  OK   ' if ok else '  FAIL ') + msg)
    if not ok:
        probs.append(msg)

def dump(ds, name):
    print(f'\n===== {name}')
    print('  dims:', dict(ds.sizes))
    for v in ds.variables:
        enc = {k: ds[v].encoding[k] for k in ('dtype', '_FillValue', 'zlib', 'complevel') if k in ds[v].encoding}
        print(f'  var {v}{ds[v].dims} {ds[v].dtype} attrs={dict(ds[v].attrs)} enc={enc}')
    print('  global attributes:')
    for k, val in ds.attrs.items():
        print(f'    {k}: {val}')

def wording(ds):
    hits = []
    for k, val in list(ds.attrs.items()) + [(f'{v}.{a}', x) for v in ds.variables for a, x in ds[v].attrs.items()]:
        s = str(val)
        if BARE.search(s.replace('called noBB in files delivered before Oct 2026', '')):
            hits.append(('bare noBB', k))
        for m in re.finditer(r'[^.;]*\bfires?\b[^.;]*', s):
            hits.append(('fire wording', f'{k}: "{m.group(0).strip()}"'))
    return hits

files = sorted(FOLDER.glob('*.nc'))
print('files:', [f.name for f in files])
U = [f for f in files if 'MDA8O3' in f.name]; Hs = [f for f in files if 'UTChourlySurfO3' in f.name]
chk(len(U) == 1 and len(Hs) == 2, f'{len(U)} MDA8 file + {len(Hs)} hourly files')

# ---- unified MDA8 ----
u = op(U[0]); dump(u, U[0].name)
print('\n--- MDA8 checks')
chk(list(map(str, u.scenario.values)) == ['BASE', 'noAnthro', 'noBBCONUS', 'noBBGlobal'], f'scenarios {list(map(str, u.scenario.values))}')
days = pd.DatetimeIndex(u.time.values)
want = pd.DatetimeIndex(np.concatenate([pd.date_range(f'{y}-04-01', f'{y}-10-31') for y in (2022, 2023)]))
chk(days.equals(want), f'{len(days)} dates {days[0].date()}..{days[-1].date()}')
chk(np.array_equal(u.lat.values, np.arange(24, 51)) and np.array_equal(u.lon.values, np.arange(-125, -65)),
    f'grid lat {u.lat.values[0]}..{u.lat.values[-1]} ({u.sizes["lat"]}), lon {u.lon.values[0]}..{u.lon.values[-1]} ({u.sizes["lon"]})')
chk(u.MDA8O3.attrs.get('units') == 'ppb', f"units {u.MDA8O3.attrs.get('units')}")
for kind, where in wording(u):
    print(f'  NOTE {kind}: {where}')

# ---- drop-in compatibility with July ----
print('\n--- compatibility with the July file', JULY.name)
j = op(JULY)
chk(list(j.data_vars) == list(u.data_vars), f'data variables July {list(j.data_vars)} / new {list(u.data_vars)}')
chk(j.MDA8O3.dims == u.MDA8O3.dims, f'dimension order {u.MDA8O3.dims}')
chk(j.MDA8O3.dtype == u.MDA8O3.dtype, f'dtype {u.MDA8O3.dtype}')
for c in ('time', 'lat', 'lon'):
    chk(np.array_equal(j[c].values, u[c].values) and j[c].dtype == u[c].dtype, f'{c} values and dtype identical ({u[c].dtype})')
chk(str(j.MDA8O3.encoding.get('_FillValue')) == str(u.MDA8O3.encoding.get('_FillValue')),
    f"_FillValue July {j.MDA8O3.encoding.get('_FillValue')} / new {u.MDA8O3.encoding.get('_FillValue')}")
print('  July scenarios:', list(map(str, j.scenario.values)), '| scenario dtype July', j.scenario.dtype, 'new', u.scenario.dtype)
# July MDA8 used DST + windows ending 07-23 + no truncation; the new one is Appendix U,
# so values differ everywhere: report how much (daily and seasonal-mean, and for the
# scenario differences Olivia uses).
v = u.MDA8O3.values; v = v[np.isfinite(v)]
chk(np.all(v == np.floor(v)), 'new MDA8 values are whole ppb (Appendix U truncation)')
chk('local standard time' in u.attrs.get('local_time', '') and 'Appendix U' in u.attrs.get('mda8_method', ''),
    'local_time / mda8_method describe Appendix U, LST')
for jn, nn in [('BASE', 'BASE'), ('noAnthro', 'noAnthro'), ('noBB', 'noBBCONUS')]:
    for yr in ('2022', '2023'):
        a = j.MDA8O3.sel(scenario=jn, time=yr).values; b = u.MDA8O3.sel(scenario=nn, time=yr).values
        d = b - a; ad = np.abs(d); sm = np.nanmean(b, 0) - np.nanmean(a, 0)
        print(f'  new {nn}{yr} - July {jn}{yr}: daily mean {np.nanmean(d):+.2f} ppb, |d|>1 ppb {100*np.nanmean(ad > 1):.0f}%, '
              f'>5 ppb {100*np.nanmean(ad > 5):.1f}%, max |d| {np.nanmax(ad):.1f}; seasonal-mean change '
              f'{np.nanmin(sm):+.2f}..{np.nanmax(sm):+.2f} (box mean {np.nanmean(sm):+.2f}) ppb; '
              f'valid-day count July {int(np.isfinite(a).sum())} / new {int(np.isfinite(b).sum())}')
for yr in ('2022', '2023'):
    dj = j.MDA8O3.sel(scenario='BASE', time=yr).mean('time') - j.MDA8O3.sel(scenario='noBB', time=yr).mean('time')
    dn = u.MDA8O3.sel(scenario='BASE', time=yr).mean('time') - u.MDA8O3.sel(scenario='noBBCONUS', time=yr).mean('time')
    print(f'  seasonal BASE-noBB(CONUS) {yr}: July box mean {float(dj.mean()):.2f}, new {float(dn.mean()):.2f} ppb; '
          f'max |change| {float(np.abs(dn - dj).max()):.2f} ppb')

# ---- hourly ----
print('\n--- hourly files')
for h in Hs:
    d = op(h); dump(d, h.name)
    yr = int(re.search(r'noBBGlobal(\d{4})', h.name).group(1))
    t = pd.DatetimeIndex(d.time.values)
    chk(list(d.data_vars) == ['O3'] and d.O3.dims == ('time', 'lat', 'lon') and d.O3.attrs.get('units') == 'ppb',
        f'{h.name}: O3(time,lat,lon) ppb')
    chk(t[0] == pd.Timestamp(f'{yr}-04-01 01:00') and t[-1] == pd.Timestamp(f'{yr}-10-31 23:00') and len(t) == 5135,
        f'time {t[0]} .. {t[-1]} UTC, {len(t)} hours')
    chk(np.array_equal(d.lat.values, u.lat.values) and np.array_equal(d.lon.values, u.lon.values), 'same grid as MDA8 file')
    chk(d.attrs.get('scenario') == 'noBBGlobal' and int(d.attrs.get('year')) == yr and
        d.attrs.get('case_name', '').endswith('noBBemisGlobalY20220401TY20231231'),
        f"scenario {d.attrs.get('scenario')} year {d.attrs.get('year')} case ...{d.attrs.get('case_name','')[-36:]}")
    for kind, where in wording(d):
        print(f'  NOTE {kind}: {where}')
print(f"\n{'ALL CHECKS PASSED' if not probs else f'{len(probs)} PROBLEM(S): ' + '; '.join(probs)}")
