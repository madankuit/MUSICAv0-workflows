#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compare_hourly_versions.py

Checks regenerated hourly gridded files (new stamp) against the earlier ones:
O3 values and time/lat/lon axes identical; reports the scenario label and case
recorded inside each file. Old 'noBB' files are matched to new 'noBBCONUS' ones.

    python compare_hourly_versions.py <new cdate> [experiment] [target]   # default ne30 1x1

    13 Oct 2026: VERSION 1.0
"""
import sys, pathlib, glob
import numpy as np, xarray as xr
_ROOT = next(p for p in pathlib.Path(__file__).resolve().parents if (p / 'config' / 'paths.py').exists())
sys.path.insert(0, str(_ROOT))
import config  # noqa: F401
from config.paths import BGO3_EXPERIMENTS
NEW = sys.argv[1]; EX = sys.argv[2] if len(sys.argv) > 2 else 'ne30'; TGT = sys.argv[3] if len(sys.argv) > 3 else '1x1'
E = BGO3_EXPERIMENTS[EX]; tag = E['targets'][TGT][2]; lp = E['label_prefix']; ok = True
for (sc, yr), cn in E['cases'].items():
    new = glob.glob(str(E['hourly_dir'] / f'{tag}_UTChourlySurfO3.{lp}{sc}{yr}.*_c{NEW}.nc'))
    old = [f for f in glob.glob(str(E['hourly_dir'] / f"{tag}_UTChourlySurfO3.{lp}{ {'noBBCONUS': 'noBB'}.get(sc, sc) if EX == 'ne30' else sc}{yr}.*.nc"))
           if not f.endswith(f'_c{NEW}.nc')]
    assert len(new) == 1 and len(old) == 1, (sc, yr, new, old)
    n = xr.open_dataset(new[0]); o = xr.open_dataset(old[0])
    same = all(np.array_equal(n[v].values, o[v].values) for v in ('time', 'lat', 'lon')) and \
        np.array_equal(n.O3.values, o.O3.values, equal_nan=True)
    print(f"{lp}{sc}{yr}: {pathlib.Path(old[0]).name} -> {pathlib.Path(new[0]).name}: values {'IDENTICAL' if same else 'DIFFERENT'}; "
          f"scenario attr {o.attrs.get('scenario')} -> {n.attrs.get('scenario')}; case ok {n.attrs['case_name'] == cn}")
    ok &= same and n.attrs['case_name'] == cn and n.attrs['scenario'] == sc
print('ALL IDENTICAL' if ok else 'PROBLEM'); sys.exit(0 if ok else 1)
