# grids/target_grids/

Regular lat/lon grids used as the **destination** of ESMF regridding from the
spectral-element grids (`ne30np4`, `ne0CONUSne30x8`).

Unlike `ne30np4/` and `ne0CONUSne30x8/`, this directory holds no workflow
scripts — only the grid descriptions that those workflows regrid *onto*.

See [`grid_files/README.md`](grid_files/README.md) for the file list and the
config constants that address them.

| Grid | Config constant | Used by |
|------|-----------------|---------|
| 1°×1° global FV | `FV_GRIDINFO_1X1` | `CONUSBGO3/regridding/` (ne30np4 → 1° CONUS) |
| 0.15°×0.15° FV | `FV_GRIDINFO_015` | `grids/ne0CONUSne30x8/regridding/` (TROPOMI comparison) |
| 0.1°×0.1° CAMS | `FV_GRIDINFO_01_CAMS` | emissions regridding |
