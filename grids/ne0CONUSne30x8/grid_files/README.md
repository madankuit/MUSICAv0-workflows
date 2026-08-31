# grid_files/

NetCDF grid description, mask, and weight files for the ne0CONUSne30x8 grid and related grids.

## ne0CONUSne30x8 Files

| File | Size | Description |
|------|------|-------------|
| `ne0CONUS_ne30x8_np4_SCRIP.nc` | ~33 MB | SCRIP grid description for ne0CONUSne30x8 (174,098 columns). Contains `grid_center_lat`, `grid_center_lon`, `grid_corner_lat/lon`, `grid_area`. Required by ESMF regridding and SE utility functions. |
| `ne0CONUSne30x8_np4_CONUSlandMaskedFalse_80kmBuffer.nc` | ~1.5 MB | Boolean mask: True = land within CONUS + 80 km coastal buffer, on ne0CONUSne30x8 ncol. |
| `US_region_masks_ne0CONUSne30x8_ncol174098.nc` | ~3.7 MB | Boolean masks for standard US sub-regions (WestCoast, Mountain, Midwest, Southwest, Southeast, Northeast) on ne0CONUSne30x8 ncol. |
| `mask_CanadaProvinces_ne0CONUS_ne30x8.nc` | ~85 MB | Province-level masks for Canadian regions on ne0CONUSne30x8 ncol. Used for biomass burning emission scaling experiments. |
| `ne0CONUSne30x8_Canadianfire_region_masks.nc` | ~8.7 MB | Regional masks for Canadian wildfire source regions (2023) on ne0CONUSne30x8 ncol. |

## ne30np4 Files

The ne30np4 SCRIP grid and its CONUS masks now live in
[`../../ne30np4/grid_files/`](../../ne30np4/grid_files/), alongside the other
ne30np4 material. Reference them as `SCRIP_NE30NP4`, `MASK_NE30NP4_CONUS_80KM`
and `MASK_NE30NP4_LOWER48_50KM` from `config/paths.py`.

## f09 Reference Files

| File | Size | Description |
|------|------|-------------|
| `f09_CONUSlandMaskedFalse_80kmBuffer.nc` | ~63 KB | CONUS land mask with 80 km buffer on the f09 (0.9°×1.25° FV) grid. Reference for finite-volume comparisons. |

---

## Usage

Take these from `config/paths.py` rather than opening them by a hard-coded path —
that way the code works from any directory and any checkout:

```python
import sys, pathlib
_ROOT = next(p for p in pathlib.Path(__file__).resolve().parents
             if (p / 'config' / 'paths.py').exists())
sys.path.insert(0, str(_ROOT))
from config.paths import SCRIP_NE0CONUSNE30X8, MASK_NE0CONUS_CONUS_80KM

import xarray as xr
scrip_ds = xr.open_dataset(SCRIP_NE0CONUSNE30X8)
mask_ds  = xr.open_dataset(MASK_NE0CONUS_CONUS_80KM)
conus_mask = mask_ds['mask']   # shape: (ncol,), dtype: bool
```

## Config constants

| Config constant | File |
|-----------------|------|
| `SCRIP_NE0CONUSNE30X8` | `ne0CONUS_ne30x8_np4_SCRIP.nc` |
| `MASK_NE0CONUS_CONUS_80KM` | `ne0CONUSne30x8_np4_CONUSlandMaskedFalse_80kmBuffer.nc` |
| `MASK_NE0CONUS_US_REGIONS` | `US_region_masks_ne0CONUSne30x8_ncol174098.nc` |
| `MASK_NE0CONUS_CANADA_PROVINCES` | `mask_CanadaProvinces_ne0CONUS_ne30x8.nc` |
| `MASK_NE0CONUS_CANADIAN_FIRE` | `ne0CONUSne30x8_Canadianfire_region_masks.nc` |
| `MASK_F09_CONUS_80KM` | `f09_CONUSlandMaskedFalse_80kmBuffer.nc` |
