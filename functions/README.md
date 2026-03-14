# functions/

General-purpose utility functions for MUSICAv0 analysis. These are **grid-agnostic** — they work with any SE grid (ne30np4, ne0CONUSne30x8, etc.) by accepting SCRIP file paths or grid metadata as inputs, or they operate entirely on arrays/datasets without grid assumptions.

Import these into grid-specific scripts rather than duplicating logic.

---

## Contents

### Emissions & Calculation

| File | Description |
|------|-------------|
| `Calc_Emis.py` | Calculate total or regional emissions from model output. Supports both FV and SE grids; accepts SCRIP file for area weighting. Returns emission totals with unit conversion (molecules/cm²/s → kg/m²/s or Tg/yr). |
| `NEIfileDatetime_shift_BYweekofday.py` | Shift NEI emission file timestamps by matching day-of-week across years. Used to map 2017 NEI files to target simulation years. |

### Regridding

| File | Description |
|------|-------------|
| `Regridding_ESMF_MTv1.py` | General ESMF regridding framework. Supports bilinear and mass-conservative (patch) methods for SE → regular lat/lon. Called by grid-specific regridding scripts. |

### SE Grid Utilities

| File | Description |
|------|-------------|
| `SE_analysis.py` | Utilities for spectral-element grids. `get_scrip_slice()` extracts a regional subset using SCRIP grid info; `get_site_index()` finds the nearest `ncol` index for a given lat/lon site. |
| `ZonalMeridional.py` | `SE_ZonalMeridional()`: area-weighted zonal or meridional means from 2D/3D SE fields using SCRIP cell areas. |

### Plotting

| File | Description |
|------|-------------|
| `Plot_2D.py` | 2D map plotting class for both FV and SE model outputs. Features: multiple projections (Cartopy), custom colormaps, log/symlog scaling, difference maps, region masking, and colorbar control. |

### Model I/O

| File | Description |
|------|-------------|
| `func_ReadMUSICAOutput.py` | Read MUSICAv0 history files. `BasicRead_ds()` for simple loading; `Read_histh1_Runds()` / `Read_histh2_Runds()` for time/level-sliced extraction from h1 (daily) and h2 (hourly) output files. |

### Model Evaluation — Statistics

| File | Description |
|------|-------------|
| `func_ModelEval_statistical_tests.py` | Statistical metrics for model–observation comparison: Pearson/Spearman correlation, Reduced Major Axis regression, and one-sided Wilcoxon signed-rank test. |

### Regional Definitions

| File | Description |
|------|-------------|
| `func_MUSICA_DefineRegion.py` | Lat/lon bounds for standard US regions (WestCoast, Mountain, Midwest, Southwest, Southeast, Northeast, CONUS). Also provides UTC timezone offsets and ~1:30 PM local time equivalents for each region. |

### Vertical Regridding (TROPOMI ↔ MUSICA)

| File | Description |
|------|-------------|
| `func_VerticalRegrid_TROPOMIAK_toMUSICAlevs.py` | Regrid TROPOMI L2 averaging kernels (TM5 pressure grid) to MUSICA hybrid pressure levels. Output on 0.1°×0.1° lat/lon. Supports NO₂ and HCHO. |
| `func_VerticalRegrid_TROPOMIAK_toMUSICAlevs_015latlon.py` | Same as above but outputs on 0.15°×0.15° lat/lon grid (matches conservative regridding output). |
| `func_VerticalRegrid_TotalCO_TROPOMIAK_toMUSICAlevs_015latlon.py` | Total CO column variant of the above. |

### Utilities

| File | Description |
|------|-------------|
| `find_missing_files.py` | Identify missing hourly files within a datetime range. `find_missing_files_v1()` for NEI-format timestamps; `find_missing_files_v2()` for CAMS-format timestamps. |

---

## Notes

- `Regridding_ESMF_MTv1.py` requires `esmpy` (the Python interface to ESMF). Install via conda: `conda install -c conda-forge esmpy`.
- `Plot_2D.py` requires `cartopy`. Install via conda: `conda install -c conda-forge cartopy`.
- All functions accept `xarray.Dataset` or `numpy` arrays; output types match inputs where possible.
