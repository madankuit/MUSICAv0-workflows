# MUSICAv0 Workflows

Analysis and postprocessing code for **MUSICAv0** (CESM2 with CAM-chem in spectral-element mode), with a focus on regional air quality studies over the continental United States.

Scripts here cover the full modeling lifecycle:

**experiment setup → emissions processing → model output postprocessing → evaluation → visualization**

---

## Repository Structure

```
MUSICAv0-workflows/
├── functions/                   # General-purpose, grid-agnostic utility functions
└── grids/
    ├── ne0CONUSne30x8/          # Variable-resolution grid (~14 km over CONUS)
    │   ├── grid_files/          # SCRIP grids, land/region masks (.nc)
    │   ├── regridding/          # Regrid SE outputs to regular lat/lon
    │   ├── emissions/           # Emissions preprocessing and sensitivity experiments
    │   ├── postprocessing/      # Model output processing (VCD, surface extraction, merging)
    │   ├── model_evaluation/    # Model–observation comparison (AQS, SLAMS, MDA8 O3)
    │   ├── satellite_comparison/# TROPOMI column comparison with averaging kernels
    │   └── plotting/            # Map visualization scripts
    └── ne30np4/                 # Standard-resolution global grid (~1°)
        └── grid_files/          # SCRIP grids and regridding weight files
```

### `functions/`
Reusable utility functions that are **grid-agnostic** or work across multiple SE grids (ne30, ne0CONUSne30x8, etc.). Import these from any analysis script. See [functions/README.md](functions/README.md).

### `grids/ne0CONUSne30x8/`
**Grid-specific** workflow scripts for the `ne0CONUSne30x8` variable-resolution CAM-SE grid (174,098 unstructured `ncol` columns; ~14 km over CONUS, ~111 km globally). Each subfolder targets a distinct stage of analysis. See [grids/ne0CONUSne30x8/README.md](grids/ne0CONUSne30x8/README.md).

### `grids/ne30np4/`
Placeholder for workflow scripts for the `ne30np4` standard-resolution global CAM-SE grid (48,602 `ncol` columns; ~111 km globally). Used for global simulations and as a lateral boundary condition source for ne0CONUS runs. See [grids/ne30np4/README.md](grids/ne30np4/README.md).

---

## Grid Overview

| Grid | Columns (`ncol`) | Resolution | Primary use |
|------|-----------------|------------|-------------|
| `ne0CONUSne30x8` | 174,098 | ~14 km CONUS / ~111 km global | Regional air quality |
| `ne30np4` | 48,602 | ~111 km global | Global / LBC source runs |
| `f09` | — | 0.9°×1.25° FV | Finite-volume reference |

SCRIP and mask files for ne0CONUSne30x8 live in [grids/ne0CONUSne30x8/grid_files/](grids/ne0CONUSne30x8/grid_files/); ne30np4 grid files in [grids/ne30np4/grid_files/](grids/ne30np4/grid_files/).

---

## Dependencies

```yaml
python >= 3.10
xarray
numpy
pandas
matplotlib
cartopy
scipy
geopandas
esmpy          # for ESMF-based conservative regridding
regionmask
```

See [environment.yml](environment.yml) for a reproducible conda environment.

---

## Data

Model output, emissions files, and observational datasets are **not** included in this repository. Scripts reference user-supplied input paths. Grid support files (SCRIP, masks) are provided as `.nc` files in `grid_files/`.

---

## Author

M. Tao
Postdoctoral Associate, MIT EAPS
Atmospheric chemistry modeling · Satellite data analysis · Regional air quality
