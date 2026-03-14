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

## ne30np4 Files (Global Background Grid)

These files are useful for lateral boundary condition (LBC) workflows and global analyses that feed into ne0CONUSne30x8 simulations.

| File | Size | Description |
|------|------|-------------|
| `ne30np4_091226_pentagons.nc` | ~5.2 MB | SCRIP grid description for the ne30np4 global grid (97,481 columns). |
| `ne30np4_091226_pentagons_CONUSlandMaskedFalse_80kmBuffer.nc` | ~434 KB | CONUS land mask with 80 km buffer on ne30np4 ncol. |
| `ne30np4_091226_pentagons_Lower48StatesCoastal50kmMaskedFalse.nc` | ~434 KB | Lower 48 states coastal mask (50 km buffer) on ne30np4 ncol. |

## f09 Reference Files

| File | Size | Description |
|------|------|-------------|
| `f09_CONUSlandMaskedFalse_80kmBuffer.nc` | ~63 KB | CONUS land mask with 80 km buffer on the f09 (0.9°×1.25° FV) grid. Reference for finite-volume comparisons. |

---

## Usage

```python
import xarray as xr

# Load SCRIP grid for se utility functions
scrip_ds = xr.open_dataset('grid_files/ne0CONUS_ne30x8_np4_SCRIP.nc')

# Load CONUS land mask
mask_ds = xr.open_dataset('grid_files/ne0CONUSne30x8_np4_CONUSlandMaskedFalse_80kmBuffer.nc')
conus_mask = mask_ds['mask']   # shape: (ncol,), dtype: bool
```
