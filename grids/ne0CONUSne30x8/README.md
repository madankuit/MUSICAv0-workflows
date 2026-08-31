# grids/ne0CONUSne30x8/

Workflow scripts for the **ne0CONUSne30x8** variable-resolution CAM-SE grid.

## Grid Specifications

| Property | Value |
|----------|-------|
| Grid type | Spectral element (SE), unstructured |
| Columns (`ncol`) | 174,098 |
| CONUS refinement | ~14 km (~ne30×8) |
| Global background | ~111 km (ne30) |
| SCRIP file | `grid_files/ne0CONUS_ne30x8_np4_SCRIP.nc` |
| Coordinate | `ncol` (1D unstructured); no lat/lon dimensions |

Output fields have shape `(time, lev, ncol)` or `(time, ncol)` — no regular lat/lon structure. Regridding is required before standard lat/lon analysis.

---

## Subfolders

### `grid_files/`
SCRIP grid description, CONUS land masks, regional masks, and Canada province masks in NetCDF format. Required by regridding scripts and SE utility functions.

### `regridding/`
Regrid ne0CONUSne30x8 SE outputs to regular lat/lon grids (0.125° linear or 0.15° conservative) for further analysis and satellite comparison.

### `emissions/`
Preprocess and modify gridded emissions for sensitivity experiments. Includes NEI-to-CAMS date alignment and regional scaling of anthropogenic or biomass burning emissions.

### `postprocessing/`
Process model history files: compute tropospheric and total VCDs (HCHO, NO₂, CO), merge hourly surface-layer outputs, and extract species at specific locations.

### `model_evaluation/`
Compare model output against surface observations (AQS/SLAMS network). Match monitor locations to grid columns, extract time series, compute MDA8 O₃, and export matched datasets to CSV/NetCDF.

### `satellite_comparison/`
Compare model columns against TROPOMI L2 retrievals using averaging kernels. Compute bias metrics, Spearman correlation, and mask NaN regions consistently between datasets. For NO₂ the stored total-column AK is converted to a tropospheric AK (`AK_trop = AK_total × AMF_total / AMF_trop`); see that subfolder's README.

### `plotting/`
Diagnostic maps for the ne0CONUSne30x8 domain: absolute and relative differences between model runs, regional overlays, and scenario comparisons.

---

## Typical Workflow

```
1. Regrid SE output to lat/lon         → regridding/
2. Compute VCD / surface extractions   → postprocessing/
3. Apply TROPOMI averaging kernels     → satellite_comparison/
4. Match with AQS observations         → model_evaluation/
5. Plot maps and diagnostics           → plotting/ + functions/Plot_2D.py
```

---

## General Functions and paths

The scripts here import from the top-level [`functions/`](../../functions/)
directory — but **do not add it to `sys.path` yourself, and never hard-code a
path to it.** Importing `config` does that for you, so shared utilities are
available by plain module name:

```python
import sys, pathlib
_ROOT = next(p for p in pathlib.Path(__file__).resolve().parents
             if (p / 'config' / 'paths.py').exists())
sys.path.insert(0, str(_ROOT))
import config                       # also puts functions/ on sys.path
from config.paths import ARCHIVE, SCRIP_NE0CONUSNE30X8

from SE_analysis import get_site_index      # from functions/
```

That bootstrap walks up to find the repo root, so it works from any
subdirectory and any checkout location. Every path these scripts use comes from
[`config/paths.py`](../../config/paths.py) — see the repository
[README](../../README.md#2-configure-your-paths).
