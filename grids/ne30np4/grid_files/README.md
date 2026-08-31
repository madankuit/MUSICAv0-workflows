# grid_files/

Grid description and mask files for the **ne30np4** spectral-element grid
(~111 km global, 48,602 `ncol` columns). These ship with the repository and are
addressed repo-relatively — no configuration needed.

---

## Files in this directory

| File | Size | Config constant | Description |
|------|------|-----------------|-------------|
| `ne30np4_091226_pentagons.nc` | ~5.2 MB | `SCRIP_NE30NP4` | SCRIP grid description for ne30np4 (48,602 columns): `grid_center_lat/lon`, `grid_corner_lat/lon`, `grid_area`. Required by ESMF regridding and the SE utility functions. |
| `ne30np4_091226_pentagons_CONUSlandMaskedFalse_80kmBuffer.nc` | ~434 KB | `MASK_NE30NP4_CONUS_80KM` | Boolean mask: CONUS land plus an 80 km buffer, on ne30np4 `ncol`. Used to zero emissions over CONUS in the CONUSBGO3 experiments. |
| `ne30np4_091226_pentagons_Lower48StatesCoastal50kmMaskedFalse.nc` | ~434 KB | `MASK_NE30NP4_LOWER48_50KM` | Lower-48 coastal mask with a 50 km buffer, on ne30np4 `ncol`. |

---

## Usage

Never open these by hard-coded path — take them from the config:

```python
import sys, pathlib
_ROOT = next(p for p in pathlib.Path(__file__).resolve().parents
             if (p / 'config' / 'paths.py').exists())
sys.path.insert(0, str(_ROOT))
from config.paths import SCRIP_NE30NP4, MASK_NE30NP4_CONUS_80KM

import xarray as xr
scrip_ds = xr.open_dataset(SCRIP_NE30NP4)
mask_ds  = xr.open_dataset(MASK_NE30NP4_CONUS_80KM)
```

---

## Related files kept outside the repository

Large derived ESMF **weight** files are not committed (3.6 MB–165 MB, and
regenerable). They resolve through `config/paths.py` to
`$DATA_ROOT/CESM22/grids/`:

| Config constant | File |
|-----------------|------|
| `WEIGHTS_NE30_TO_1X1` | `ESMFmap_ne30np4_TO_1x1_conserve_c20260708.nc` |
| `WEIGHTS_FV09X125_TO_NE30` | `ESMFmap_0.9x1.25_ne30np4_cubit_conserve_cams.nc` |
| `SCRIP_NE30NP4_ALT` | `ne30np4_grid_c20241105.nc` (alternative ne30np4 SCRIP) |

Regenerate the ne30→1° weights with
[`CONUSBGO3/regridding/gen_ne30_to_1x1_weights.py`](../../../CONUSBGO3/regridding/gen_ne30_to_1x1_weights.py).
