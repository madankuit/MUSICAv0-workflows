# CONUSBGO3/regridding/archive/

Method-development scripts from the ne30np4 → 1° regrid work (8 Jul 2026),
kept for provenance. **Neither produces the deliverable** — the production
pipeline is [`../Regrid_ne30_surfO3_to_1x1_conserve.py`](../Regrid_ne30_surfO3_to_1x1_conserve.py).

Both operate on a single trial day (2022-07-20) for BASE / noAnthro / noBB and
cross-check the regridded field against the raw ne30 column values at a sample
of monitors, so the two are directly comparable.

| Script | Method | Outcome |
|--------|--------|---------|
| `apply_conserve_trial.py` | ESMF first-order **conservative** remap, applied as a sparse mat-mul of the offline weights | **Adopted.** Became `../Regrid_ne30_surfO3_to_1x1_conserve.py`. |
| `trial_regrid_linear.py` | **Linear (Delaunay) interpolation** between SE column centres via `LinearNDInterpolator` | **Rejected.** Does not conserve mass. |

## Why conservative won

The deliverable is meant to be area-integrated and compared against other
gridded products, so the remap has to conserve the field's mass. Linear
interpolation between spectral-element column centres does not: it is smooth
and cheap, but the cell means it produces do not integrate back to the source
field. First-order conservative remapping does, and its destination weight-sums
verify to 1 (checked in `../gen_ne30_to_1x1_weights.py`).

The two scripts also differ in target grid — the linear trial used an ad-hoc
CONUS box (`lon 210–310`, `lat 0–70`), whereas the conservative path uses the
global 1° FV grid description shipped at
`grids/target_grids/grid_files/FV1x1grid_info_c20241105.nc`, which carries the
cell **corners** that conservative remapping requires.

## Still useful

`apply_conserve_trial.py` doubles as a quick smoke test: if the ne30→1° weight
file is ever regenerated, run it to confirm the operator and the monitor
cross-check still behave before launching the full multi-year job.

Both read every path from `config/paths.py` like the rest of the repo.
