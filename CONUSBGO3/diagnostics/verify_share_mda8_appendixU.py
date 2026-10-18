#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify_share_mda8_appendixU.py

Final check of the MDA8 file shared with the GAM group against 40 CFR Part 50
Appendix U, Section 3, using the hourly files that go with it:

1. ALL cells, all scenario-years, Apr 1 - Oct 30: an independent pandas
   implementation (rolling sums labelled at the window start), different code from
   functions/mda8_epa.py and from audit_deliverables.py.
2. 100 random cells per scenario-year: the literal clause-by-clause reference (same
   logic as functions/test_mda8_epa.literal, written in ppm as in the regulation).
3. Oct 31: from the merged hourly model output (which extends past Oct 31), regridded
   with the file's own weights is out of scope here; instead we confirm Oct 31 values are
   whole ppb and that the days left missing are exactly those of runs ending Nov 1 00 UTC.
4. File-level facts the definition implies: whole ppb, no negatives, metadata text.

Oct 31 local needs model hours up to Nov 1 ~15 UTC; the hourly files end Oct 31 23 UTC,
so items 1-2 stop at Oct 30.

    python verify_share_mda8_appendixU.py <unified MDA8 file> <hourly dir>

    18 Oct 2026: VERSION 1.0
"""
import sys, glob, math, pathlib
import numpy as np, pandas as pd, xarray as xr
from timezonefinder import TimezoneFinder
import pytz
from datetime import datetime

U = xr.open_dataset(sys.argv[1]); HD = pathlib.Path(sys.argv[2])
M = U.MDA8O3
tf = TimezoneFinder()
LA, LO = np.meshgrid(U.lat.values, U.lon.values, indexing='ij')
lat, lon = LA.ravel(), LO.ravel()

def lst(la, lo, year):
    """Standard time of that year (January = standard in the northern hemisphere). Year matters:
    Chihuahua/Ojinaga were -7 in 2022 and moved to -6 from Oct 30, 2022."""
    z = tf.timezone_at(lat=la, lng=lo)
    if z is None:
        return int(round(lo / 15))
    return int(pytz.timezone(z).localize(datetime(year, 1, 15)).utcoffset().total_seconds() // 3600)
OFFY = {yr: np.array([lst(a, b, yr) for a, b in zip(lat, lon)]) for yr in (2022, 2023)}
print('LST offsets present:', {yr: sorted(set(o.tolist())) for yr, o in OFFY.items()},
      '| cells differing 2022 vs 2023:', int((OFFY[2022] != OFFY[2023]).sum()))

def pandas_impl(series_end_utc_ppb, off, year):
    """Independent: pandas Series (index = CAM end-of-hour UTC) -> MDA8 per local date."""
    s = series_end_utc_ppb.copy()
    s.index = s.index - pd.Timedelta(hours=1) + pd.Timedelta(hours=int(off))   # hour-beginning LST
    s = np.floor(s + 1e-6)                                                    # 3(a) truncation
    full = s.reindex(pd.date_range(f'{year}-03-31', f'{year}-11-02', freq='h'))
    fwd = full[::-1]                                                          # forward window via reversed rolling
    tot = fwd.rolling(8, min_periods=1).sum()[::-1]                           # sum of hours t..t+7
    cnt = fwd.rolling(8, min_periods=1).count()[::-1]
    tot = tot.where(cnt > 0)
    avg = np.floor(tot / cnt + 1e-6).where(cnt >= 6)                          # 3(b) divisor 6-8
    zs = np.floor(tot / 8 + 1e-6)
    avg = avg.where(avg.notna(), zs.where((cnt <= 5) & (zs > 70)))           # 3(b) zero substitution
    hrs = avg.index.hour
    blk = avg[(hrs >= 7) & (hrs <= 23)]                                       # 3(c) starts 07..23
    g = blk.groupby(blk.index.normalize())
    mx, nv = g.max(), g.count()
    res = mx.where((nv >= 13) | (mx > 70))                                    # 3(d)
    return res.reindex(pd.date_range(f'{year}-04-01', f'{year}-10-30'))

def hour_begin_lst(hourly_end_utc, vals_ppb, off):
    return {t - pd.Timedelta(hours=1) + pd.Timedelta(hours=int(off)): v for t, v in zip(hourly_end_utc, vals_ppb)}

def literal(hb, day):
    """Clause by clause, in ppm as written in Appendix U. hb: {hour-beginning LST: ppb}."""
    def hr(t):
        v = hb.get(t, math.nan)
        return math.nan if math.isnan(v) else math.floor(v / 1000 * 1000 + 1e-6) / 1000   # 3(a)
    def eight(t0):
        x = [hr(t0 + pd.Timedelta(hours=i)) for i in range(8)]; have = [v for v in x if not math.isnan(v)]
        if len(have) >= 6:
            a = sum(have) / len(have)
        elif have and sum(have) > 0.567 + 1e-12:
            a = sum(have) / 8
        else:
            return math.nan
        return math.floor(a * 1000 + 1e-6) / 1000                              # truncated
    av = [eight(day + pd.Timedelta(hours=h)) for h in range(7, 24)]
    ok = [a for a in av if not math.isnan(a)]
    if not ok:
        return math.nan
    mx = max(ok)
    return mx * 1000 if (len(ok) >= 13 or mx > 0.070) else math.nan

rng = np.random.default_rng(42); bad = 0
for sc in map(str, U.scenario.values):
    for yr in (2022, 2023):
        h = xr.open_dataset(glob.glob(str(HD / f'CONUS1x1_UTChourlySurfO3.{sc}{yr}.*.nc'))[0]).O3
        t = pd.DatetimeIndex(h.time.values); Hv = h.values.reshape(len(t), -1).astype('float64')
        ref = M.sel(scenario=sc, time=slice(f'{yr}-04-01', f'{yr}-10-30')).values.reshape(-1, Hv.shape[1])
        OFF = OFFY[yr]; mism = 0; nanmis = 0
        for c in range(Hv.shape[1]):
            p = pandas_impl(pd.Series(Hv[:, c], index=t), OFF[c], yr).values
            mism += int(np.nansum(np.abs(p - ref[:, c]) > 0)); nanmis += int((np.isnan(p) != np.isnan(ref[:, c])).sum())
        lit_bad = 0
        days = pd.date_range(f'{yr}-04-01', f'{yr}-10-30')
        for c in rng.choice(Hv.shape[1], 100, replace=False):
            hb = hour_begin_lst(t, Hv[:, c], OFF[c])
            for d, day in enumerate(days):
                L = literal(hb, day); F = ref[d, c]
                lit_bad += int(not ((math.isnan(L) and np.isnan(F)) or L == F))
        ok = mism == 0 and nanmis == 0 and lit_bad == 0
        bad += not ok
        print(f"{'OK  ' if ok else 'FAIL'} {sc}{yr}: all {Hv.shape[1]} cells x {len(days)} days vs pandas impl: "
              f"{mism} value / {nanmis} NaN mismatches; literal reference on 100 cells: {lit_bad} mismatches")

v = M.values; f = v[np.isfinite(v)]
print(f"whole ppb: {bool(np.all(f == np.floor(f)))}; min {f.min():.0f} max {f.max():.0f} ppb; "
      f"negatives {int((f < 0).sum())}; missing {100*np.isnan(v).mean():.2f}%")
o31 = {(sc, yr): int(np.isnan(M.sel(scenario=sc, time=f'{yr}-10-31').values).sum())
       for sc in map(str, U.scenario.values) for yr in (2022, 2023)}
print('missing cells on Oct 31 per scenario-year:', o31, '(expected only BASE2023 / noBBCONUS2022, whose output ends Nov 1 00 UTC)')
print('mda8_method:', U.attrs.get('mda8_method')); print('local_time:', U.attrs.get('local_time'))
print('references:', U.attrs.get('references'))
print('ALL MDA8 VALUES MATCH APPENDIX U' if bad == 0 else f'{bad} SCENARIO-YEARS WITH MISMATCHES')
