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

## Getting Started

### 1. Clone and set up the environment

```bash
git clone https://github.com/madankuit/MUSICAv0-workflows.git
cd MUSICAv0-workflows
conda env create -f environment.yml
conda activate musica-workflows
```

### 2. Configure your paths

All cluster paths (model output directories, grid files, regridding weights, etc.) are collected in one place:

```
svante_MUSICA_paths.py
```

This file contains **one person's** Svante paths as a reference. Copy it and update the paths to match your own directories:

```bash
cp svante_MUSICA_paths.py my_paths.py   # or keep as svante_MUSICA_paths.py locally
```

Then edit the variables at the top of the file (paths to CESM archive, regridding weights, output save directories, etc.).

To verify all your paths exist on the cluster, run:

```bash
python check_svante_paths.py            # checks all entries
python check_svante_paths.py --fix      # writes svante_MUSICA_paths_updated.py with status annotations
```

### 3. Set the `functions/` path in each script

Scripts in `grids/` import shared utilities from `functions/`. Each script has a **USER CONFIGURATION** block near the top — set `repo_root` (or `functions_path`) there:

```python
# --- USER CONFIGURATION ---
repo_root = '/path/to/MUSICAv0-workflows'
import sys; sys.path.insert(0, repo_root + '/functions')
```

### 4. Typical analysis workflow (ne0CONUSne30x8)

```
1. Regrid SE output to lat/lon         →  grids/ne0CONUSne30x8/regridding/
2. Compute VCDs / surface extractions  →  grids/ne0CONUSne30x8/postprocessing/
3. Apply TROPOMI averaging kernels     →  grids/ne0CONUSne30x8/satellite_comparison/
4. Match with AQS/SLAMS observations   →  grids/ne0CONUSne30x8/model_evaluation/
5. Plot maps and diagnostics           →  grids/ne0CONUSne30x8/plotting/
```

Each subfolder has its own README with script-level descriptions.

> **Note (TROPOMI NO₂ averaging kernel).** In step 3, the TROPOMI NO₂ averaging kernel
> stored in the L2 product is the *total-column* kernel and is converted to a tropospheric
> kernel via `AK_trop = AK_total × (AMF_total / AMF_trop)` (TROPOMI ATBD S5P-KNMI-L2-0005-RP),
> matching the L2 match/recalc convention. This requires regridded `AMF_total`/`AMF_trop`
> fields; see `functions/README.md` and `grids/ne0CONUSne30x8/satellite_comparison/README.md`.

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
