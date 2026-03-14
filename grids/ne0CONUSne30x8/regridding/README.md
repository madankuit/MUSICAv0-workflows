# regridding/

Regrid ne0CONUSne30x8 spectral-element (SE) output from unstructured `ncol` coordinates to regular lat/lon grids. Two approaches are provided:

1. **Linear interpolation** — fast, suitable for diagnostic maps
2. **Mass-conservative (ESMF)** — required for column-integrated quantities and satellite comparisons

---

## Files

### Linear Regridding

| File | Output resolution | Description |
|------|------------------|-------------|
| `Regrid_SE_015FV.py` | 1°×1° (ne30) or 0.125°×0.125° (ne0CONUSne30x8) | Linear interpolation of SE surface-layer output to regular lat/lon using `scipy.interpolate`. Contains grid-specific functions for ne30 and ne0CONUSne30x8. Outputs NetCDF files. |

### Conservative Regridding (ESMF, for TROPOMI overpass)

| File | Output resolution | Description |
|------|------------------|-------------|
| `ConservativeRegrid_TROPOMITime_MUSICAOutputs.py` | 0.15°×0.15° | Mass-conservative ESMF regrid of hourly ne0CONUSne30x8 SE outputs at times near the TROPOMI overpass (~1:30 PM LT). Processes vertical and non-vertical variables separately. |
| `ConservativeRegrid_MergeAllOneVerticalVars.py` | — | Merge regridded non-vertical variables (e.g., `TROP_P`, `PS`) across dates after running `ConservativeRegrid_TROPOMITime_MUSICAOutputs.py`. |
| `ConservativeRegrid_MergeAllLev_eachDate.py` | — | Merge regridded vertical variables (all levels) for each date after conservative regridding. |

---

## Workflow

```
1. Run ConservativeRegrid_TROPOMITime_MUSICAOutputs.py
   → outputs individual date/variable files at 0.15°×0.15°

2. Run ConservativeRegrid_MergeAllLev_eachDate.py
   → assembles vertical variables by date

3. Run ConservativeRegrid_MergeAllOneVerticalVars.py
   → assembles non-vertical variables across all dates

4. Use merged files in satellite_comparison/ for VCD + TROPOMI AK application
```

---

## Notes

- Conservative regridding requires `esmpy` and the SCRIP file `grid_files/ne0CONUS_ne30x8_np4_SCRIP.nc`.
- The 0.15°×0.15° output grid matches the TM5-MP pressure coordinate used in TROPOMI L2 products.
- Linear regridding (`Regrid_SE_015FV.py`) is much faster but does not conserve mass — do not use for VCD calculations.
