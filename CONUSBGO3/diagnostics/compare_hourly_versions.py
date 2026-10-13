#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compare_hourly_versions.py

Checks regenerated ne30 hourly gridded files (new stamp) against the earlier ones:
O3 values and time/lat/lon axes identical; reports the scenario label and case
recorded inside each file. Old 'noBB' files are matched to new 'noBBCONUS' ones.

    python compare_hourly_versions.py <new cdate>

    13 Oct 2026: VERSION 1.0
"""
import sys, pathlib, glob
import numpy as np, xarray as xr
_ROOT = next(p for p in pathlib.Path(__file__).resolve().parents if (p / 'config' / 'paths.py').exists())
sys.path.insert(0, str(_ROOT))
import config  # noqa: F401
from config.paths import BGO3_EXPERIMENTS
E = BGO3_EXPERIMENTS['ne30']; NEW = sys.argv[1]; ok = True
for (sc, yr), cn in E['cases'].items():
    new = glob.glob(str(E['hourly_dir'] / f'CONUS1x1_UTChourlySurfO3.{sc}{yr}.*_c{NEW}.nc'))
    old = [f for f in glob.glob(str(E['hourly_dir'] / f"CONUS1x1_UTChourlySurfO3.{ {'noBBCONUS': 'noBB'}.get(sc, sc) }{yr}.*.nc"))
           if not f.endswith(f'_c{NEW}.nc')]
    assert len(new) == 1 and len(old) == 1, (sc, yr, new, old)
    n = xr.open_dataset(new[0]); o = xr.open_dataset(old[0])
    same = all(np.array_equal(n[v].values, o[v].values) for v in ('time', 'lat', 'lon')) and \
        np.array_equal(n.O3.values, o.O3.values, equal_nan=True)
    print(f"{sc}{yr}: {pathlib.Path(old[0]).name} -> {pathlib.Path(new[0]).name}: values {'IDENTICAL' if same else 'DIFFERENT'}; "
          f"scenario attr {o.attrs.get('scenario')} -> {n.attrs.get('scenario')}; case ok {n.attrs['case_name'] == cn}")
    ok &= same and n.attrs['case_name'] == cn and n.attrs['scenario'] == sc
print('ALL IDENTICAL' if ok else 'PROBLEM'); sys.exit(0 if ok else 1)
