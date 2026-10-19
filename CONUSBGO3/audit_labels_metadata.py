#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
audit_labels_metadata.py

Label / metadata audit of the CONUSBGO3 share set (Oct 2026). For every current
deliverable it checks that the file name, the labels inside the file and the
simulation it came from all agree, and that no outdated wording is left:

- name: scenario + year label (ne0CONUSne30x8_ prefix on that grid), grid tag,
  version stamp c20261012 for gridded files
- attributes: scenario, year, case_name (against a hand-written list of the
  simulations), grid, model description (inventory + met of that setup),
  regrid weight / source SCRIP / destination grid matching the grid and target
- unified MDA8: scenario coordinate, case_names entries, scenario text, title,
  local-time rule, time axis
- no bare 'noBB' label anywhere (only the documented "called noBB before" note)
- units and dimensions (803 monitors; ppb vs mol/mol)

It also lists every other .nc in the output folders, marking files the share set
supersedes. Prints one line per file; exits 1 if anything is wrong.

    python audit_labels_metadata.py [--list-files]   # --list-files: print the share set paths only

    18 Oct 2026: VERSION 1.1 - MDA8 stamp c20261018 (Appendix U), hourly c20261012
    13 Oct 2026: VERSION 1.0
"""
import sys, pathlib, glob, re
import numpy as np, pandas as pd, xarray as xr

_ROOT = next(p for p in pathlib.Path(__file__).resolve().parents if (p / 'config' / 'paths.py').exists())
sys.path.insert(0, str(_ROOT))
import config  # noqa: F401
from config.paths import BGO3_EXPERIMENTS
# --experiments=ne30[,ne0CONUSne30x8]: audit only these sets (ne0 outputs deleted 2026-10-18)
_ex=[a.split('=',1)[1].split(',') for a in sys.argv if a.startswith('--experiments=')]
EXPS={k:v for k,v in BGO3_EXPERIMENTS.items() if not _ex or k in _ex[0]}
sys.argv=[a for a in sys.argv if not a.startswith('--experiments=')]


STAMP = '20261018'          # MDA8 (unified + point MDA8 files): Appendix U recompute
HSTAMP = '20261018'         # hourly gridded files (re-issued 2026-10-19 with the MDA8 stamp; values unchanged)
SIM = {   # written out by hand
    'ne30': {('BASE', 2022): 'f.e22.FCnudged.ne30_ne30_mg17.BGO3.BASEY20220401TY20230401',
             ('BASE', 2023): 'f.e22.FCnudged.ne30_ne30_mg17.BGO3.BASEY20230401TY20231101',
             ('noAnthro', 2022): 'f.e22.FCnudged.ne30_ne30_mg17.BGO3.noANTHROemisCONUS80kmBufferY20220401TY20221101',
             ('noAnthro', 2023): 'f.e22.FCnudged.ne30_ne30_mg17.BGO3.noANTHROemisCONUS80kmBufferY20230401TY20231101',
             ('noBBCONUS', 2022): 'f.e22.FCnudged.ne30_ne30_mg17.BGO3.noBBemisCONUS80kmBufferY20220401TY20221101',
             ('noBBCONUS', 2023): 'f.e22.FCnudged.ne30_ne30_mg17.BGO3.noBBemisCONUS80kmBufferY20230401TY20231101',
             ('noBBGlobal', 2022): 'f.e22.FCnudged.ne30_ne30_mg17.BGO3.noBBemisGlobalY20220401TY20231231',
             ('noBBGlobal', 2023): 'f.e22.FCnudged.ne30_ne30_mg17.BGO3.noBBemisGlobalY20220401TY20231231'},
    'ne0CONUSne30x8': {
             ('BASE', 2023): 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.hourlyNEIemisCONUS.From20230101.01',
             ('BASE', 2024): 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.hourlyNEIemisCONUS.From20230101.01',
             ('noBBCONUS', 2023): 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.hourlyNEInoBBemisCONUS.From20230101.02',
             ('noBBCONUS', 2024): 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.hourlyNEInoBBemisCONUS.From20230101.02'},
}
GRIDNAME = {'ne30': 'ne30np4', 'ne0CONUSne30x8': 'ne0CONUSne30x8'}
MODELKEY = {'ne30': ['ne30np4', 'MERRA-2', 'CAMS-GLOB-ANT v6.2'],
            'ne0CONUSne30x8': ['ne0CONUSne30x8', 'GEOS-5', 'NEI2022v2', 'NOT the same setup']}
BARE_NOBB = re.compile(r'noBB(?![A-Za-z])')
ALLOWED_NOBB = 'called noBB in files delivered before Oct 2026'
problems = []; share = []


def report(path, issues, info=''):
    tag = 'OK  ' if not issues else 'FAIL'
    print(f"{tag} {path.name}" + (f"  [{info}]" if info else '') + ''.join(f"\n       - {i}" for i in issues))
    if issues:
        problems.append((path.name, issues))


def text_issues(ds):
    """Bare 'noBB' anywhere in global/variable attributes or the scenario coordinate."""
    out = []
    blobs = {f'attr {k}': str(v) for k, v in ds.attrs.items()}
    for v in ds.variables:
        blobs.update({f'{v}.{k}': str(x) for k, x in ds[v].attrs.items()})
    if 'scenario' in ds.coords:
        blobs['scenario coord'] = ' '.join(map(str, ds.scenario.values))
    for k, s in blobs.items():
        if BARE_NOBB.search(s.replace(ALLOWED_NOBB, '')):
            out.append(f"bare 'noBB' in {k}: ...{s[max(0, BARE_NOBB.search(s).start()-40):BARE_NOBB.search(s).start()+40]}...")
    return out


def common(ds, ex, sc, yr, kind):
    iss = []
    if ds.attrs.get('case_name') != SIM[ex][(sc, yr)]:
        iss.append(f"case_name {ds.attrs.get('case_name')} != {SIM[ex][(sc, yr)]}")
    if str(ds.attrs.get('scenario')) != sc:
        iss.append(f"scenario attr {ds.attrs.get('scenario')} != {sc}")
    if int(ds.attrs.get('year', -1)) != yr:
        iss.append(f"year attr {ds.attrs.get('year')} != {yr}")
    src = ds.attrs.get('source', ds.attrs.get('model', ''))
    miss = [k for k in MODELKEY[ex] if k not in src and not (ex == 'ne30' and k == 'ne30np4' and 'ne30np4' in src)]
    if kind == 'points':
        if ex == 'ne0CONUSne30x8' and ds.attrs.get('grid') != 'ne0CONUSne30x8':
            iss.append(f"grid attr {ds.attrs.get('grid')}")
        if ex == 'ne30' and ds.attrs.get('grid') != 'ne30np4':
            iss.append(f"grid attr {ds.attrs.get('grid')}")
    if miss:
        iss.append(f"model/source text lacks {miss}")
    mf = str(ds.attrs.get('merged_input_file', ''))
    if not mf.startswith(SIM[ex][(sc, yr)] + '.cam.h2.surflev.O3.' + f'{yr}-04-01'):
        iss.append(f"merged_input_file {mf}")
    return iss


for ex, E in EXPS.items():
    lp = E['label_prefix']
    # ---- unified MDA8 ----
    for tgt, (gi, wgt, tag) in E['targets'].items():
        f = pathlib.Path(E['regrid_dir']) / (E['unified_stem'].format(tag=tag) + f'_c{STAMP}.nc')
        share.append(f)
        if not f.exists():
            report(f, ['missing']); continue
        ds = xr.open_dataset(f); iss = text_issues(ds)
        want = [s for s in E['scenarios']]
        if list(map(str, ds.scenario.values)) != want:
            iss.append(f"scenario coord {list(ds.scenario.values)}")
        for (sc, yr), cn in SIM[ex].items():
            if f"{lp}{sc}{yr}: {cn}" not in ds.attrs.get('case_names', ''):
                iss.append(f"case_names lacks {lp}{sc}{yr}: {cn}")
        res = '1x1 deg' if tgt == '1x1' else '0.15x0.15 deg'
        for k, need in [('title', GRIDNAME[ex]), ('title', res), ('horizontal_regrid', res),
                        ('regrid_weight_file', pathlib.Path(wgt).name), ('regrid_dest_grid', pathlib.Path(gi).name),
                        ('regrid_source_scrip', pathlib.Path(E['scrip']).name), ('local_time', 'local standard time'),
                        ('mda8_method', 'Appendix U'), ('mda8_method', '07:00-23:00 LST'), ('mda8_values', 'whole ppb'),
                        ('references', '40 CFR Part 50, Appendix U'), ('references', '80 FR 65292'),
                        ('scenario_dim', ', '.join(want)), ('scenarios', 'noBBCONUS')]:
            if need not in str(ds.attrs.get(k, '')):
                iss.append(f"{k} lacks '{need}'")
        for k in MODELKEY[ex]:
            if k not in ds.attrs.get('source', ''):
                iss.append(f"source lacks '{k}'")
        days = pd.DatetimeIndex(ds.time.values)
        exp_days = pd.DatetimeIndex(np.concatenate([pd.date_range(f'{y}-04-01', f'{y}-10-31') for y in E['years']]))
        if not days.equals(exp_days):
            iss.append('time axis')
        if ds.MDA8O3.attrs.get('units') != 'ppb':
            iss.append('MDA8O3 units')
        _v = ds.MDA8O3.values; _v = _v[np.isfinite(_v)]
        if not np.all(_v == np.floor(_v)):
            iss.append('MDA8 not whole ppb (Appendix U truncation)')
        report(f, iss, f"{ex} {tgt}: scenarios {want}, {len(days)} days")
    # ---- hourly gridded ----
    for tgt, (gi, wgt, tag) in E['targets'].items():
        for (sc, yr), cn in SIM[ex].items():
            hits = sorted(glob.glob(str(E['hourly_dir'] / f'{tag}_UTChourlySurfO3.{lp}{sc}{yr}.*.nc')))
            cur = [h for h in hits if h.endswith(f'_c{HSTAMP}.nc')]
            if not cur:
                cur = hits[-1:]                         # ne0 hourly were made c20261010 (unchanged by later fixes)
            f = pathlib.Path(cur[0]) if cur else pathlib.Path(f'{tag}_UTChourlySurfO3.{lp}{sc}{yr}.MISSING')
            share.append(f)
            if not cur:
                report(f, ['missing']); continue
            ds = xr.open_dataset(f); iss = text_issues(ds) + common(ds, ex, sc, yr, 'hourly')
            res = '1x1 deg' if tgt == '1x1' else '0.15x0.15 deg'
            for need in [GRIDNAME[ex], res, f'{lp}{sc}{yr}']:
                if need not in ds.attrs.get('title', ''):
                    iss.append(f"title lacks '{need}'")
            if pathlib.Path(wgt).name not in ds.attrs.get('regrid_weight_file', ''):
                iss.append('regrid_weight_file')
            t = pd.DatetimeIndex(ds.time.values)
            if t[0] != pd.Timestamp(f'{yr}-04-01 01:00') or t[-1] != pd.Timestamp(f'{yr}-10-31 23:00'):
                iss.append(f'time {t[0]}..{t[-1]}')
            if ds.O3.attrs.get('units') != 'ppb':
                iss.append('O3 units')
            stamp = re.search(r'_c(\d{8})\.nc$', f.name).group(1)
            report(f, iss, f"{ex} {tgt} {sc}{yr}, stamp c{stamp}")
    # ---- points ----
    ptag = 'All803' if ex == 'ne30' else ''
    for (sc, yr), cn in SIM[ex].items():
        for kind, var, unit in [('LocalTimeMDA8O3', 'MDA8O3', 'ppb'), ('UTChourlyO3', 'O3', 'mol/mol')]:
            hits = glob.glob(str(E['points_dir'] / f'{kind}.{lp}{sc}{yr}.GivenMonitors{ptag}.*.nc'))
            if kind == 'LocalTimeMDA8O3':
                hits = [h for h in hits if h.endswith(f'_c{STAMP}.nc')]
            f = pathlib.Path(hits[0]) if len(hits) == 1 else pathlib.Path(f'{kind}.{lp}{sc}{yr}.GivenMonitors{ptag}.MISSING')
            share.append(f)
            if len(hits) != 1:
                report(f, [f'{len(hits)} matches']); continue
            ds = xr.open_dataset(f); iss = text_issues(ds) + common(ds, ex, sc, yr, 'points')
            if ds.sizes.get('AQS_code') != 803:
                iss.append(f"AQS_code {ds.sizes.get('AQS_code')}")
            if ds[var].attrs.get('units') != unit:
                iss.append(f"{var} units {ds[var].attrs.get('units')}")
            if kind == 'LocalTimeMDA8O3':
                if 'local standard time' not in ds.attrs.get('mda8_method', ''):
                    iss.append('mda8_method lacks local standard time')
                if '40 CFR Part 50, Appendix U' not in ds.attrs.get('references', ''):
                    iss.append('references attribute missing')
                v = ds[var].values; v = v[np.isfinite(v)]
                if not np.all(v == np.floor(v)):
                    iss.append('MDA8 not whole ppb')
            ci = ds.attrs.get('column_index_file', '')
            want_ci = 'MatchedMonitors_ne30_ColIdx_exact.csv' if ex == 'ne30' else 'MatchedMonitors_ne0CONUSne30x8_ColIdx.csv'
            if ci != want_ci:
                iss.append(f"column_index_file {ci} (want {want_ci})")
            report(f, iss, f"{ex} {sc}{yr} {kind}")

if '--list-files' in sys.argv:
    print('\n'.join(str(p) for p in share)); sys.exit(0)

# ---- everything else in the output folders ----
print('\n=== other .nc files in the output folders (not in the share set) ===')
shareset = {str(p) for p in share}
for ex, E in EXPS.items():
    for d in sorted({E['regrid_dir'], E['hourly_dir'], E['points_dir']}):
        for f in sorted(glob.glob(str(pathlib.Path(d) / '*.nc'))):
            if f not in shareset:
                n = pathlib.Path(f).name
                why = ('old noBB label' if BARE_NOBB.search(n) else
                       'legacy monitor matching (781 sites)' if '.GivenMonitors.' in n and ex == 'ne30' else
                       'pre-Appendix-U MDA8 / older stamp' if 'MDA8O3' in n else
                       'older stamp, values identical to c20261012' if 'UTChourlySurfO3' in n else '?')
                print(f"   superseded: {d.name if hasattr(d,'name') else d}/{n}  <- {why}")
print(f"\n{'ALL LABELS AND METADATA CONSISTENT' if not problems else f'{len(problems)} FILE(S) WITH PROBLEMS'}")
sys.exit(1 if problems else 0)
