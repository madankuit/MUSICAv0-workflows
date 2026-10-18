#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
add_mda8_references.py

Attribute-only edit (Oct 2026): adds the 'references' global attribute (EPA
Appendix U and the 2015 O3 NAAQS rule, functions/mda8_epa.REFERENCES) to the c<stamp>
MDA8 files (unified gridded + point) made before the attribute was in the code.
Data untouched: md5 of every variable's values checked before/after.

    python add_mda8_references.py 20261018 [--dry-run]

    18 Oct 2026: VERSION 1.0
"""
import sys, pathlib, hashlib
import numpy as np, netCDF4
_ROOT = next(p for p in pathlib.Path(__file__).resolve().parents if (p / 'config' / 'paths.py').exists())
sys.path.insert(0, str(_ROOT))
import config  # noqa: F401
from config.paths import BGO3_EXPERIMENTS
from mda8_epa import REFERENCES

STAMP = sys.argv[1]; DRY = '--dry-run' in sys.argv

def data_md5(nc):
    h = hashlib.md5()
    for v in sorted(nc.variables):
        a = np.asarray(nc[v][:]); h.update(v.encode())
        h.update('\x00'.join(map(str, a.ravel())).encode() if a.dtype.kind in 'OUS' else
                 np.ma.filled(nc[v][:], np.nan).tobytes() if np.ma.isMaskedArray(nc[v][:]) else a.tobytes())
    return h.hexdigest()

files = []
for E in BGO3_EXPERIMENTS.values():
    files += sorted(pathlib.Path(E['regrid_dir']).glob(f'MUSICAv0_*_MDA8O3_*_c{STAMP}.nc'))
    files += sorted(pathlib.Path(E['points_dir']).glob(f'LocalTimeMDA8O3.*_c{STAMP}.nc'))
bad = 0
for f in files:
    with netCDF4.Dataset(f) as nc:
        before = data_md5(nc); had = nc.getncattr('references') if 'references' in nc.ncattrs() else None
    if not DRY and had != REFERENCES:
        with netCDF4.Dataset(f, 'a') as nc:
            nc.setncattr('references', REFERENCES)
    with netCDF4.Dataset(f) as nc:
        after = data_md5(nc); now = nc.getncattr('references') if 'references' in nc.ncattrs() else None
    ok = before == after and (DRY or now == REFERENCES)
    bad += not ok
    print(f"{'OK  ' if ok else 'FAIL'} {f.name}: references {'already set' if had == REFERENCES else 'added' if not DRY else 'to add'}; "
          f"data md5 {'unchanged' if before == after else 'CHANGED'}")
print(f'{len(files)} files, {bad} problems{" (dry run)" if DRY else ""}')
sys.exit(1 if bad or not files else 0)
