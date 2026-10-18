#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
audit_deliverables.py

Independent audit of the CONUSBGO3 deliverables (Oct 2026). Nothing here reuses
the pipeline's merge / matching / regrid / MDA8 code: values are recomputed from
the RAW h2 history files with separate implementations and compared with the
delivered files.

A. Inventory: every deliverable's dims, dates, scenarios, units, ranges, NaN%,
   and that its case_name / case_names attributes name the right simulation.
B. Points: for sample monitors (incl. corrected and newly added ones) the column
   is found with shapely point-in-polygon on the SCRIP grid; hourly O3 is read
   from the raw h2 files; MDA8 is recomputed with plain numpy. Compared with the
   point files for every (scenario, year), ne30 (All803) and ne0CONUSne30x8.
C. Regrid: for sample target cells, weights are rebuilt from shapely polygon
   intersections (lon-lat plane, cos(lat) area factor) instead of the ESMF file;
   raw h2 O3 at sample hours is averaged and compared with the hourly gridded files.
D. Grid MDA8: MDA8 recomputed (numpy) from the hourly gridded files at sample cells
   and compared with the unified MDA8 files (Apr 1 - Oct 30; the hourly files stop
   at Oct 31 23 UTC, so Oct 31 local is not complete in them).

MDA8 convention checked: 40 CFR 50 Appendix U, written independently here (loops)
from functions/mda8_epa.py: local standard time, hour-beginning = CAM stamp - 1 h,
hourly and 8-h values truncated to whole ppb, 17 periods starting 07-23 LST, >=6 h
(or zero-substituted > 70 ppb), day valid with >=13 or maximum > 70 ppb.

    python audit_deliverables.py [cdate YYYYMMDD]

    18 Oct 2026: VERSION 1.1 - Appendix U MDA8 (independent implementation)
    12 Oct 2026: VERSION 1.0
"""
import sys, pathlib, glob, re
import numpy as np, pandas as pd, xarray as xr
import netCDF4
from shapely.geometry import Point, Polygon, box
from timezonefinder import TimezoneFinder
import pytz
from datetime import datetime

_ROOT = next(p for p in pathlib.Path(__file__).resolve().parents
             if (p / 'config' / 'paths.py').exists())
sys.path.insert(0, str(_ROOT))
import config  # noqa: F401
from config.paths import BGO3_EXPERIMENTS, case_hist_dir
# --experiments=ne30[,ne0CONUSne30x8]: audit only these sets (ne0 outputs deleted 2026-10-18)
_ex=[a.split('=',1)[1].split(',') for a in sys.argv if a.startswith('--experiments=')]
EXPS={k:v for k,v in BGO3_EXPERIMENTS.items() if not _ex or k in _ex[0]}
sys.argv=[a for a in sys.argv if not a.startswith('--experiments=')]


RNG = np.random.default_rng(20261012)
CDATE = sys.argv[1] if len(sys.argv) > 1 else '20261018'   # MDA8 (unified + point) stamp to audit
FAIL = []
def check(ok, msg):
    print(("   OK   " if ok else "   FAIL ") + msg)
    if not ok:
        FAIL.append(msg)

# expected simulations, written out by hand (not taken from the pipeline)
EXPECT = {
    'ne30': {'BASE': {2022: 'BASEY20220401TY20230401', 2023: 'BASEY20230401TY20231101'},
             'noAnthro': {2022: 'noANTHROemisCONUS80kmBufferY20220401TY20221101',
                          2023: 'noANTHROemisCONUS80kmBufferY20230401TY20231101'},
             'noBBCONUS': {2022: 'noBBemisCONUS80kmBufferY20220401TY20221101',
                           2023: 'noBBemisCONUS80kmBufferY20230401TY20231101'},
             'noBBGlobal': {2022: 'noBBemisGlobalY20220401TY20231231', 2023: 'noBBemisGlobalY20220401TY20231231'}},
    'ne0CONUSne30x8': {'BASE': {2023: 'hourlyNEIemisCONUS.From20230101.01', 2024: 'hourlyNEIemisCONUS.From20230101.01'},
                 'noBBCONUS': {2023: 'hourlyNEInoBBemisCONUS.From20230101.02',
                               2024: 'hourlyNEInoBBemisCONUS.From20230101.02'}},
}
for ex, d in EXPECT.items():
    for sc, yy in d.items():
        for yr, short in yy.items():
            assert BGO3_EXPERIMENTS[ex]['cases'][(sc, yr)].endswith(short), (ex, sc, yr)

tf = TimezoneFinder()
def summer_offset(lat, lon, year):
    """Local STANDARD time offset (h), written independently of functions/mda8_epa.py.
    (Name kept for the call sites; Appendix U uses LST.)"""
    tz = tf.timezone_at(lat=lat, lng=lon)
    if tz is None:
        return int(round(lon / 15.0))
    jan = pytz.timezone(tz).localize(datetime(year, 1, 15))          # northern-hemisphere winter = standard time
    return int(jan.utcoffset().total_seconds() // 3600)

def mda8_numpy(times_utc, vals, off, year):
    """Independent Appendix-U MDA8 (plain loops): Series local-date -> whole-ppb value or NaN.
    times_utc are CAM end-of-hour stamps -> hour-beginning LST = stamp - 1 h + off."""
    hb = pd.DatetimeIndex(times_utc) - pd.Timedelta(hours=1) + pd.Timedelta(hours=off)
    h = pd.Series(np.floor(np.asarray(vals, dtype='float64') + 1e-6), index=hb)   # hourly truncation
    h = h[~h.index.duplicated()]
    out = []
    days = pd.date_range(f'{year}-04-01', f'{year}-10-31')
    for d in days:
        wins = []
        for st in range(7, 24):                                         # periods starting 07..23 LST
            t0 = d + pd.Timedelta(hours=st)
            x = h.reindex(pd.date_range(t0, t0 + pd.Timedelta(hours=7), freq='h')).values
            n = int(np.isfinite(x).sum())
            if n >= 6:
                a = np.floor(np.nansum(x) / n + 1e-6)
            elif n > 0 and np.floor(np.nansum(x) / 8 + 1e-6) > 70:          # zeros substituted
                a = np.floor(np.nansum(x) / 8 + 1e-6)
            else:
                a = np.nan
            wins.append(a)
        wins = np.array(wins); nv = int(np.isfinite(wins).sum())
        mx = np.nanmax(wins) if nv else np.nan
        out.append(mx if (nv >= 13 or (nv and mx > 70)) else np.nan)
    return pd.Series(out, index=days)

def raw_files(case, year, end_mmdd='11-01'):
    rx = re.compile(r'\.cam\.h2\.(\d{4}-\d{2}-\d{2})-\d{5}\.nc$')
    fs = []
    for f in glob.glob(str(case_hist_dir(case) / f'{case}.cam.h2.*.nc')):
        d = pd.Timestamp(rx.search(f).group(1))
        if pd.Timestamp(f'{year}-04-01') <= d <= pd.Timestamp(f'{year}-{end_mmdd}'):
            fs.append((d, f))
    return [f for _, f in sorted(fs)]

def raw_series(case, year, cols):
    """Hourly O3 (mol/mol) at the lowest level for columns `cols`, from raw h2."""
    T, V = [], []
    for f in raw_files(case, year):
        with netCDF4.Dataset(f) as nc:
            t = netCDF4.num2date(nc['time'][:], nc['time'].units, nc['time'].calendar)   # cftime (noleap)
            T.extend(pd.to_datetime([str(x) for x in t]))
            V.append(np.asarray(nc['O3'][:, -1, :])[:, cols])
    T = pd.DatetimeIndex(T); V = np.concatenate(V)
    keep = T >= pd.Timestamp(f'{year}-04-01 01:00')            # common first stamp
    T, V = T[keep], V[keep]
    _, first = np.unique(T, return_index=True)
    return T[first], V[first]

def contains_col(scrip, lat, lon):
    clat = scrip['grid_center_lat'].values; clon = scrip['grid_center_lon'].values
    d = np.abs(clat - lat) + np.abs(((clon - lon + 180) % 360) - 180)
    for i in np.argsort(d)[:60]:
        xs = lon + (((scrip['grid_corner_lon'].values[i] - lon + 180) % 360) - 180)
        if Polygon(list(zip(xs, scrip['grid_corner_lat'].values[i]))).buffer(0).contains(Point(lon, lat)):
            return int(i)
    return None

# ======================================================================= A
print("\n=== A. inventory ===")
inv = []
for ex, E in EXPS.items():
    lp = E['label_prefix']
    for tgt, (_, _, tag) in E['targets'].items():
        stem = E['unified_stem'].format(tag=tag)
        f = sorted(glob.glob(str(E['regrid_dir'] / f'{stem}_c{CDATE}.nc')))
        check(len(f) == 1, f'{ex} {tgt}: unified file {stem}_c{CDATE}.nc present')
        if not f:
            continue
        m = xr.open_dataset(f[0]); M = m.MDA8O3
        check(list(m.scenario.values) == list(EXPECT[ex]), f'{ex} {tgt}: scenarios {list(m.scenario.values)}')
        days = pd.DatetimeIndex(m.time.values)
        want = pd.DatetimeIndex(np.concatenate([pd.date_range(f'{y}-04-01', f'{y}-10-31') for y in E['years']]))
        check(days.equals(want), f'{ex} {tgt}: 428 dates Apr1-Oct31 of {E["years"]}')
        # units sanity only (mol/mol or a double 1e9 would be off by orders of magnitude); the model
        # has real fire-plume extremes up to ~460 ppb at 0.15 deg, reported, not failed
        check(M.attrs.get('units') == 'ppb' and 0 <= float(M.min()) and float(M.max()) < 1000,   # 0 = truncated <1 ppb (App. U)
              f'{ex} {tgt}: MDA8 range {float(M.min()):.1f}..{float(M.max()):.1f} ppb, NaN {100*float(np.isnan(M).mean()):.2f}%, '
              f'cell-days >150 ppb {100*float((M > 150).mean()):.3f}%')
        for sc, yy in EXPECT[ex].items():
            for yr, short in yy.items():
                check(f'{lp}{sc}{yr}: ' in m.attrs['case_names'] and short in
                      m.attrs['case_names'].split(f'{lp}{sc}{yr}: ')[1].split(';')[0],
                      f'{ex} {tgt}: case_names {sc}{yr} -> ...{short}')
        inv.append((ex, tgt, f[0]))
    tagpt = 'All803' if ex == 'ne30' else ''
    for sc, yy in EXPECT[ex].items():
        for yr, short in yy.items():
            for kind in ['LocalTimeMDA8O3', 'UTChourlyO3']:
                f = glob.glob(str(E['points_dir'] / f'{kind}.{lp}{sc}{yr}.GivenMonitors{tagpt}.*.nc'))
                ok = len(f) == 1
                if ok:
                    p = xr.open_dataset(f[0]); v = p[list(p.data_vars)[0]]
                    ok = (p.attrs['case_name'].endswith(short) and p.sizes['AQS_code'] == 803)
                    rng = f'{float(v.min()):.3g}..{float(v.max()):.3g} {v.attrs.get("units")}'
                check(ok, f'{ex} points {kind} {lp}{sc}{yr}: 803 sites, case ...{short}' + (f', {rng}' if len(f) == 1 else ''))
    for tgt, (_, _, tag) in E['targets'].items():
        hf = sorted(glob.glob(str(E['hourly_dir'] / f'{tag}_UTChourlySurfO3.*.nc')))
        pre = f'{tag}_UTChourlySurfO3.'
        lab_of = lambda h: pathlib.Path(h).name[len(pre):].split('.')[0]      # tag may contain dots
        for sc, yy in EXPECT[ex].items():
            for yr, short in yy.items():
                lab = f'{lp}{sc}{yr}'
                alt = lab.replace('noBBCONUS', 'noBB')      # ne30 hourly files from Jul 2026 use the old label
                hit = [h for h in hf if lab_of(h) in (lab, alt)]
                ok = len(hit) == 1
                if ok:
                    ds = xr.open_dataset(hit[0]); ok = ds.attrs['case_name'].endswith(short)
                    tt = pd.DatetimeIndex(ds.time.values)
                    ok = ok and tt[0] <= pd.Timestamp(f'{yr}-04-01 01:00') and tt[-1] == pd.Timestamp(f'{yr}-10-31 23:00')
                check(ok, f'{ex} {tgt}: hourly gridded {lab}: 1 file, case ...{short}, Apr 1 .. Oct 31 23 UTC')

# ======================================================================= B
print("\n=== B. points recomputed from raw h2 (independent column lookup + numpy MDA8) ===")
SITES = {'ne30': ['aqs_340130003', 'aqs_040134010', 'aqs_060811001', 'aqs_530330080', 'aqs_230010014', 'aqs_410050004'],
         'ne0CONUSne30x8': ['aqs_340130003', 'aqs_060374009', 'aqs_530330080', 'aqs_120712002']}
for ex, E in EXPS.items():
    scrip = xr.open_dataset(E['scrip'])
    mon = pd.read_csv(E.get('colidx_exact', E['colidx']), dtype={'AQS_code': str}).set_index('AQS_code')
    extra = list(RNG.choice([a for a in mon.index if 25 < mon.lat[a] < 49 and a not in SITES[ex]], 3, replace=False))
    sites = SITES[ex] + extra
    cols = []
    for a in sites:
        c = contains_col(scrip, mon.lat[a], mon.lon[a])
        check(c == int(mon.MUSICA0_colIndex[a]), f'{ex} {a}: shapely column {c} == pipeline column {mon.MUSICA0_colIndex[a]}')
        cols.append(c)

    lp = E['label_prefix']; tagpt = 'All803' if ex == 'ne30' else ''
    for sc, yy in EXPECT[ex].items():
        for yr in yy:
            case = E['cases'][(sc, yr)]
            offs = [summer_offset(mon.lat[a], mon.lon[a], yr) for a in sites]
            T, V = raw_series(case, yr, cols)
            pm = xr.open_dataset(glob.glob(str(E['points_dir'] / f'LocalTimeMDA8O3.{lp}{sc}{yr}.GivenMonitors{tagpt}.*.nc'))[0])
            ph = xr.open_dataset(glob.glob(str(E['points_dir'] / f'UTChourlyO3.{lp}{sc}{yr}.GivenMonitors{tagpt}.*.nc'))[0])
            dmax_h = dmax_m = 0.0; nan_mis = 0; edge = 0
            # runs whose output ends before Nov 1 ~08 UTC cannot complete Oct 31 local in the West;
            # the pipeline and this audit treat the missing tail differently -> compared separately
            short_tail = T[-1] < pd.Timestamp(f'{yr}-11-01 08:00')
            for k, a in enumerate(sites):
                h = ph.O3.sel(AQS_code=a).to_series()
                com = h.index.intersection(T)
                dmax_h = max(dmax_h, float(np.nanmax(np.abs(h.reindex(com).values -
                                                            pd.Series(V[:, k], index=T).reindex(com).values))) * 1e9)
                mine = mda8_numpy(T, V[:, k] * 1e9, offs[k], yr)
                theirs = pm.MDA8O3.sel(AQS_code=a).to_series().reindex(mine.index)
                if short_tail:
                    last = mine.index == pd.Timestamp(f'{yr}-10-31')
                    edge += int(((np.abs(mine.values - theirs.values) > 1e-3) |
                                 (np.isnan(mine.values) != np.isnan(theirs.values)))[last].sum())
                    mine, theirs = mine[~last], theirs[~last]
                nan_mis += int((np.isnan(mine.values) != np.isnan(theirs.values)).sum())
                dmax_m = max(dmax_m, float(np.nanmax(np.abs(mine.values - theirs.values))))
            check(dmax_h == 0 and dmax_m < 1e-3 and nan_mis == 0,
                  f'{ex} {lp}{sc}{yr}: {len(sites)} sites, raw-vs-file hourly max|d| {dmax_h:.1e} ppb, '
                  f'MDA8 max|d| {dmax_m:.1e} ppb, NaN mismatches {nan_mis}' +
                  (f' (Apr 1-Oct 30; output ends {T[-1]:%m-%d %H} UTC, Oct 31 differs at {edge}/{len(sites)} sites '
                   f'= incomplete-record edge)' if short_tail else ''))

# ======================================================================= C
print("\n=== C. regrid recomputed with shapely overlap weights (raw h2, sample hours) ===")
for ex, E in EXPS.items():
    scrip = xr.open_dataset(E['scrip'])
    clat = scrip['grid_center_lat'].values; clon = scrip['grid_center_lon'].values
    clon180 = np.where(clon > 180, clon - 360, clon)
    for tgt, (_, _, tag) in E['targets'].items():
        res = 1.0 if tgt == '1x1' else 0.15
        lp = E['label_prefix']
        sc = 'BASE'; yr = E['years'][-1]
        hf = glob.glob(str(E['hourly_dir'] / f'{tag}_UTChourlySurfO3.{lp}{sc}{yr}.*.nc'))
        H = xr.open_dataset(sorted(hf)[-1]).O3
        pts = [(40.7, -74.2), (34.0, -118.2), (45.6, -122.6), (39.7, -105.0), (29.7, -95.1), (44.0, -114.0)]
        hours = pd.DatetimeIndex([f'{yr}-07-15 18:00', f'{yr}-08-20 06:00'])
        files = {d: [f for f in raw_files(E['cases'][(sc, yr)], yr) if f'.h2.{d:%Y-%m-%d}-' in f] for d in hours.normalize()}
        dmax = 0.0
        for la, lo in pts:
            cy = float(H.lat.sel(lat=la, method='nearest')); cx = float(H.lon.sel(lon=lo, method='nearest'))
            cell = box(cx - res / 2, cy - res / 2, cx + res / 2, cy + res / 2)
            near = np.where((np.abs(clat - cy) < 2.5) & (np.abs(clon180 - cx) < 2.5 / np.cos(np.deg2rad(cy))))[0]
            w = {}
            for i in near:
                xs = cx + (((scrip['grid_corner_lon'].values[i] - cx + 180) % 360) - 180)
                ys = scrip['grid_corner_lat'].values[i]
                a = Polygon(list(zip(xs, ys))).buffer(0).intersection(cell)
                if a.area > 0:
                    w[i] = a.area * np.cos(np.deg2rad(a.centroid.y))
            idx = np.array(list(w)); wt = np.array([w[i] for i in idx]); wt /= wt.sum()
            for t in hours:
                f = [g for g in files[t.normalize()]]
                with netCDF4.Dataset(f[0]) as nc:
                    tt = pd.to_datetime([str(x) for x in netCDF4.num2date(nc['time'][:], nc['time'].units,
                                                                          nc['time'].calendar)])
                    k = int(np.where(tt == t)[0][0])
                    v = np.asarray(nc['O3'][k, -1, :])[idx]
                mine = float((wt * v).sum()) * 1e9
                theirs = float(H.sel(time=t, lat=cy, lon=cx))
                dmax = max(dmax, abs(mine - theirs) / theirs)
        check(dmax < 0.01, f'{ex} {tgt} {lp}{sc}{yr}: {len(pts)} cells x {len(hours)} hours, '
                           f'shapely-weight vs file max rel diff {100*dmax:.3f}% (planar vs great-circle; <1% expected)')

# ======================================================================= D
print("\n=== D. unified MDA8 vs numpy MDA8 from the hourly gridded files ===")
for ex, tgt, uf in inv:
    E = BGO3_EXPERIMENTS[ex]; tag = E['targets'][tgt][2]; lp = E['label_prefix']
    U = xr.open_dataset(uf).MDA8O3
    for sc in EXPECT[ex]:
        for yr in E['years']:
            lab = f'{lp}{sc}{yr}'
            hf = glob.glob(str(E['hourly_dir'] / f'{tag}_UTChourlySurfO3.{lab}.*.nc')) or \
                 glob.glob(str(E['hourly_dir'] / f'{tag}_UTChourlySurfO3.{lab.replace("noBBCONUS", "noBB")}.*.nc'))
            H = xr.open_dataset(sorted(hf)[-1]).O3
            dmax = 0.0; nm = 0
            for _ in range(4):
                iy = int(RNG.integers(0, H.sizes['lat'])); ix = int(RNG.integers(0, H.sizes['lon']))
                la, lo = float(H.lat[iy]), float(H.lon[ix])
                s = H.isel(lat=iy, lon=ix).to_series()
                if s.isna().all():
                    continue
                mine = mda8_numpy(s.index, s.values, summer_offset(la, lo, yr), yr).loc[:f'{yr}-10-30']
                theirs = U.sel(scenario=sc, time=slice(f'{yr}-04-01', f'{yr}-10-30')).isel(lat=iy, lon=ix).to_series()
                theirs.index = mine.index
                nm += int((np.isnan(mine.values) != np.isnan(theirs.values)).sum())
                dmax = max(dmax, float(np.nanmax(np.abs(mine.values - theirs.values))))
            check(dmax < 1e-3 and nm == 0, f'{ex} {tgt} {lab}: 4 random cells, MDA8 max|d| {dmax:.1e} ppb, NaN mismatches {nm}')

print(f"\n{'ALL CHECKS PASSED' if not FAIL else f'{len(FAIL)} CHECK(S) FAILED:'}")
for f in FAIL:
    print("  -", f)
sys.exit(1 if FAIL else 0)
