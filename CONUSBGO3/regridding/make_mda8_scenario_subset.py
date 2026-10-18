#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_mda8_scenario_subset.py

Writes a unified MDA8 file holding only some scenarios of an existing one (values
copied unchanged), with the scenario description, case names and title rewritten
to cover only those scenarios. Used for the Oct 2026 share with the GAM group:
BASE + noBBGlobal only (fire O3 = BASE - noBBGlobal).

    python make_mda8_scenario_subset.py <source unified file> BASE noBBGlobal

Output: next to the source, named ..._MDA8O3_<scen1>-<scen2>_<period>_c<stamp>.nc
(same c<stamp> as the source: same data version).

    18 Oct 2026: VERSION 1.0
"""
import sys, pathlib, re
from datetime import datetime
import numpy as np, xarray as xr

DESC = {
    'BASE': 'BASE = all emissions',
    'noAnthro': 'noAnthro = CONUS anthropogenic emissions zeroed (land, 80 km buffer)',
    'noBBCONUS': 'noBBCONUS = QFED2.6 biomass-burning emissions zeroed over CONUS land + 80 km buffer only',
    'noBBGlobal': ('noBBGlobal = biomass-burning emissions zeroed globally (QFED2.6 plus CMIP6 fire DMS and '
                   'num_so4_a1), so BASE - noBBGlobal = O3 from all fires incl. transported smoke (e.g. Canadian) '
                   'and the ~1 ppb hemispheric fire background'),
}

src = pathlib.Path(sys.argv[1]); keep = sys.argv[2:]
m = xr.open_dataset(src)
assert set(keep) <= set(map(str, m.scenario.values)), (keep, m.scenario.values)
out = m.sel(scenario=keep)
assert np.array_equal(out.MDA8O3.values, m.MDA8O3.sel(scenario=keep).values, equal_nan=True)

a = dict(m.attrs)
labs = [c.strip() for c in a['case_names'].split(';')]
a['case_names'] = '; '.join(c for c in labs if re.match(r'(\w+?)(\d{4}):', c).group(1) in keep)
a['merged_input_files'] = '; '.join(c.strip() for c in a.get('merged_input_files', '').split(';')
                                    if c.strip() and re.match(r'(\w+?)(\d{4}):', c.strip()).group(1) in keep)
a['scenarios'] = '; '.join(DESC[s] for s in keep)
a['scenario_dim'] = ', '.join(keep)
a['title'] = re.sub(r'\d+ scenarios', f'{len(keep)} scenarios ({", ".join(keep)})', a['title'])
if 'noBBGlobal' in keep:
    a['scenario_notes'] = ('noBBGlobal is one continuous run from 2022-04-01, so its 2023 season starts from a '
                           'no-fire state; BASE2023 branches from BASE2022 at 2023-04-01. Grid, MERRA-2 nudging, '
                           'all other emissions and the 2022 initial condition are identical.')
a['scenario_source'] = (f"subset ({', '.join(keep)}) of {src.name}, values unchanged; "
                        f"created {datetime.now():%Y-%m-%d}")
out.attrs = a
stamp = re.search(r'_c(\d{8})\.nc$', src.name).group(1)
name = src.name.replace('_MDA8O3_', f"_MDA8O3_{'-'.join(keep)}_").replace(f'_c{stamp}.nc', f'_c{stamp}.nc')
dst = src.parent / name
assert not dst.exists(), dst
out.to_netcdf(dst, format='NETCDF4', encoding={'MDA8O3': {'zlib': True, 'complevel': 4,
                                                          '_FillValue': np.float32(np.nan)}})
chk = xr.open_dataset(dst)
assert list(map(str, chk.scenario.values)) == keep
assert np.array_equal(chk.MDA8O3.values, m.MDA8O3.sel(scenario=keep).values, equal_nan=True)
print('wrote', dst)
print('scenario:', list(chk.scenario.values), '| dims', dict(chk.MDA8O3.sizes))
for k in ['title', 'scenarios', 'scenario_dim', 'case_names', 'scenario_source']:
    print(f'  {k}: {chk.attrs[k]}')
others = [s for s in ('noBBCONUS', 'noAnthro') if s not in keep]
leak = [k for k, v in chk.attrs.items() if any(o in str(v) for o in others)]
print('  attributes still mentioning', others, ':', leak or 'none')
