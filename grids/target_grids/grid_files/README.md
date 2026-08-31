# grid_files/

Regular lat/lon **destination grids** used as the target of ESMF regridding from
the spectral-element grids. These are small grid *descriptions* (12 KB–118 KB),
so they ship with the repository and need no configuration.

---

## Files in this directory

| File | Size | Config constant | Description |
|------|------|-----------------|-------------|
| `FV1x1grid_info_c20241105.nc` | ~12 KB | `FV_GRIDINFO_1X1` | 1°×1° global finite-volume grid info with CF cell bounds. Destination for the ne30np4 → 1° conservative remap used by CONUSBGO3. |
| `FV_gridinfo_0.15_c20231204.nc` | ~71 KB | `FV_GRIDINFO_015` | 0.15°×0.15° FV grid info. Destination for the ne0CONUSne30x8 → 0.15° conservative remap used for the TROPOMI comparison. |
| `FV_gridinfo_CAMS_c20210219.nc` | ~117 KB | `FV_GRIDINFO_01_CAMS` | 0.1°×0.1° CAMS emissions grid info. Destination for emissions regridding. |

These carry cell **corners/bounds**, which first-order conservative remapping
requires — a plain lat/lon axis pair is not sufficient.

---

## Usage

```python
import sys, pathlib
_ROOT = next(p for p in pathlib.Path(__file__).resolve().parents
             if (p / 'config' / 'paths.py').exists())
sys.path.insert(0, str(_ROOT))
from config.paths import FV_GRIDINFO_1X1, SCRIP_NE30NP4

import esmpy
dstgrid = esmpy.Grid(filename=str(FV_GRIDINFO_1X1),
                     filetype=esmpy.FileFormat.GRIDSPEC,
                     add_corner_stagger=True)
```

---

## Weight files are *not* here

The ESMF **weight** files produced from these grids are large (3.6 MB–165 MB)
and derived, so they stay on the cluster and resolve through `config/paths.py`
(`WEIGHTS_NE30_TO_1X1`, `WEIGHTS_NE0CONUS_TO_015`, `WEIGHTS_NE0CONUS_TO_01`,
`WEIGHTS_FV09X125_TO_NE30`).

Conservative maps are **not reversible** — a weight file built FV→SE cannot be
used SE→FV. Generate the direction you need; see
[`CONUSBGO3/regridding/gen_ne30_to_1x1_weights.py`](../../../CONUSBGO3/regridding/gen_ne30_to_1x1_weights.py)
for a worked example (needs `esmpy`).
