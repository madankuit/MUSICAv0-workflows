#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify_noBBGlobal_addition.py

Checks the 4-scenario unified MDA8 file (Oct 2026) before it is shared:

1. BASE / noAnthro / noBBCONUS are identical to the delivered 3-scenario file
   (whose 'noBB' is noBBCONUS).
2. Every scenario's case_name attribute and the point files point at the right
   simulation (noBBCONUS -> ...noBBemisCONUS80kmBuffer..., noBBGlobal -> ...noBBemisGlobal...).
3. noBBGlobal2022 ~= BASE2022 on the first day (both branch from the same state
   on 2022-04-01) and diverges later.
4. Gridded vs point MDA8 for noBBGlobal (nearest 1x1 cell): r and bias.
5. BASE - noBBGlobal >= BASE - noBBCONUS on the CONUS mean, and the Quebec-smoke
   day 2023-06-07 at NYC.

    python verify_noBBGlobal_addition.py <delivered 3-scenario unified file>

    9 Oct 2026: VERSION 1.0
"""
import sys, pathlib, glob
import numpy as np, xarray as xr

_ROOT = next(p for p in pathlib.Path(__file__).resolve().parents
             if (p / 'config' / 'paths.py').exists())
sys.path.insert(0, str(_ROOT))
import config  # noqa: F401
from config.paths import (BGO3_CASES, BGO3_YEARS, BGO3_GIVEN_MONITORS_DIR,
                          bgo3_unified_mda8_latest)

OLD = sys.argv[1]
NEW = bgo3_unified_mda8_latest()
n = xr.open_dataset(NEW); o = xr.open_dataset(OLD)
print("new:", NEW, "\nold:", OLD)
print("scenarios:", list(n.scenario.values), " dims:", dict(n.MDA8O3.sizes))
assert list(n.scenario.values) == ["BASE", "noAnthro", "noBBCONUS", "noBBGlobal"]

# 1. copied scenarios unchanged
for new_s, old_s in [("BASE", "BASE"), ("noAnthro", "noAnthro"), ("noBBCONUS", "noBB")]:
    a = n.MDA8O3.sel(scenario=new_s).values; b = o.MDA8O3.sel(scenario=old_s).values
    same = np.array_equal(np.isnan(a), np.isnan(b)) and np.array_equal(a[~np.isnan(a)], b[~np.isnan(b)])
    print(f"1. {new_s:10s} == old {old_s:9s}: {'IDENTICAL' if same else 'DIFFERENT'}")
    assert same
assert np.array_equal(n.time.values, o.time.values) and np.array_equal(n.lat.values, o.lat.values) \
    and np.array_equal(n.lon.values, o.lon.values)

# 2. which simulation each scenario came from
print("2. case_names:", n.attrs["case_names"])
for (sc, yr), cn in BGO3_CASES.items():
    assert f"{sc}{yr}: {cn}" in n.attrs["case_names"]
assert "noBBemisCONUS80kmBuffer" in BGO3_CASES[("noBBCONUS", 2022)]
assert "noBBemisCONUS80kmBuffer" in BGO3_CASES[("noBBCONUS", 2023)]
assert all("noBBemisGlobal" in BGO3_CASES[("noBBGlobal", y)] for y in BGO3_YEARS)
print("   merged_input_files:", n.attrs["merged_input_files"])
print("   scenario_source:", n.attrs["scenario_source"])
pts = {}
for yr in BGO3_YEARS:
    f = glob.glob(str(BGO3_GIVEN_MONITORS_DIR / f"LocalTimeMDA8O3.noBBGlobal{yr}.GivenMonitors.*.nc"))
    assert len(f) == 1, f
    pts[yr] = xr.open_dataset(f[0])
    print(f"   point file {yr}: case_name={pts[yr].attrs['case_name']}  merged={pts[yr].attrs['merged_input_file']}")
    assert pts[yr].attrs["case_name"] == BGO3_CASES[("noBBGlobal", yr)]

# 3. same start, diverging later
M = n.MDA8O3
for day in ["2022-04-01", "2022-04-03", "2022-07-15"]:
    d = (M.sel(scenario="BASE", time=day) - M.sel(scenario="noBBGlobal", time=day))
    print(f"3. {day}: BASE - noBBGlobal  mean {float(d.mean()):+.3f}  max|d| {float(np.abs(d).max()):.3f} ppb")

# 4. gridded vs point (nearest cell)
for yr in BGO3_YEARS:
    p = pts[yr]
    g = M.sel(scenario="noBBGlobal", time=str(yr)).sel(
        lat=xr.DataArray(p.lat.values, dims="AQS_code"),
        lon=xr.DataArray(p.lon.values, dims="AQS_code"), method="nearest")
    pv = p.MDA8O3.transpose("AQS_code", "time").values
    gv = g.transpose("AQS_code", "time").values
    ok = np.isfinite(pv) & np.isfinite(gv)
    r = np.corrcoef(pv[ok], gv[ok])[0, 1]
    print(f"4. noBBGlobal{yr} gridded vs point: r={r:.3f}  bias(grid-point)={np.mean(gv[ok]-pv[ok]):+.2f} ppb  n={ok.sum()}")

# 5. all-fire >= US-fire attribution
base = M.sel(scenario="BASE")
dC = (base - M.sel(scenario="noBBCONUS")).mean("time"); dG = (base - M.sel(scenario="noBBGlobal")).mean("time")
print(f"5. 2-yr CONUS mean: BASE-noBBCONUS {float(dC.mean()):.2f} ppb, BASE-noBBGlobal {float(dG.mean()):.2f} ppb, "
      f"cells with dG < dC: {100*float((dG < dC).mean()):.1f}%")
nyc = M.sel(time="2023-06-07").sel(lat=41, lon=-74, method="nearest")
print("   NYC 2023-06-07 MDA8:", {str(s): round(float(v), 1) for s, v in zip(n.scenario.values, nyc.values)})
print("all checks passed")
