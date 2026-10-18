#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
relabel_bb_wording.py

Attribute-only edit (Oct 2026): the scenario descriptions in the gridded share files
said "fire(s)" where the experiments remove biomass-burning (BB) emissions (QFED2.6,
the inventory's own term). Rewrites the global attributes `scenarios` and
`scenario_notes` in place with the same wording the pipeline now writes. Data are
not touched: an md5 of every variable's values is compared before and after.

Files: every current (c20261012) unified MDA8 and hourly gridded file of all
experiments, i.e. the gridded part of audit_labels_metadata.py's share set.

    python relabel_bb_wording.py [--dry-run]

    18 Oct 2026: VERSION 1.0
"""
import sys, pathlib, re, hashlib
import numpy as np
import netCDF4
_ROOT = next(p for p in pathlib.Path(__file__).resolve().parents if (p / 'config' / 'paths.py').exists())
sys.path.insert(0, str(_ROOT))
import config  # noqa: F401
from config.paths import BGO3_EXPERIMENTS

REPL = [  # (old, new) applied to the joined attribute text, in order
    ('(fires in Canada/Alaska/Mexico kept, so BASE - noBBCONUS = O3 from US fires',
     '(biomass-burning emissions in Canada/Alaska/Mexico kept, so BASE - noBBCONUS = O3 from US biomass burning'),
    ('CMIP6 fire DMS', 'CMIP6 biomass-burning DMS'),
    ('O3 from all fires incl. transported smoke and the hemispheric fire background',
     'O3 from all biomass burning incl. transported smoke and the hemispheric biomass-burning background'),
    ('no-fire state', 'no-biomass-burning state'),
    ('so CONUS fires have been off in noBBCONUS', 'so CONUS biomass-burning emissions have been off in noBBCONUS'),
    ('non-fire emissions', 'non-biomass-burning emissions'),
]
FIRE = re.compile(r'\bfires?\b|no-fire|non-fire')
DRY = '--dry-run' in sys.argv

def data_md5(nc):
    h = hashlib.md5()
    for v in sorted(nc.variables):
        a = nc[v][:]
        h.update(v.encode())
        if np.asarray(a).dtype.kind in 'OUS':                      # text (e.g. scenario names): hash the values
            h.update('\x00'.join(map(str, np.asarray(a).ravel())).encode())
        else:
            h.update(np.ma.filled(a, np.nan).tobytes() if np.ma.isMaskedArray(a) else np.asarray(a).tobytes())
    return h.hexdigest()

files = []
for ex, E in BGO3_EXPERIMENTS.items():
    files += sorted(f for f in pathlib.Path(E['regrid_dir']).glob('*_c20261012.nc')
                    if f.name.startswith(E['unified_stem'].split('{')[0]) and '_MDA8O3_2' in f.name)   # unified only
    files += sorted(pathlib.Path(E['hourly_dir']).glob('*_c20261012.nc'))
bad = 0
for f in files:
    with netCDF4.Dataset(f, 'r') as nc:
        before = data_md5(nc); attrs = {k: nc.getncattr(k) for k in ('scenarios', 'scenario_notes') if k in nc.ncattrs()}
    new = {}
    for k, s in attrs.items():
        t = s
        for o, n in REPL:
            t = t.replace(o, n)
        new[k] = t
    left = {k: FIRE.findall(v) for k, v in new.items() if FIRE.search(v)}
    changed = [k for k in attrs if attrs[k] != new[k]]
    if not DRY and changed:
        with netCDF4.Dataset(f, 'a') as nc:
            for k in changed:
                nc.setncattr(k, new[k])
    with netCDF4.Dataset(f, 'r') as nc:
        after = data_md5(nc)
        now = {k: nc.getncattr(k) for k in attrs}
    ok = before == after and not left and (DRY or all(now[k] == new[k] for k in attrs))
    bad += not ok
    print(f"{'OK  ' if ok else 'FAIL'} {f.name}: attrs changed {changed or 'none'}; data md5 "
          f"{'unchanged' if before == after else 'CHANGED'}; remaining fire wording {left or 'none'}")
print(f"{len(files)} files, {bad} problems{' (dry run)' if DRY else ''}")
sys.exit(1 if bad else 0)
