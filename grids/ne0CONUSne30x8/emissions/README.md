# emissions/

Scripts for preprocessing and modifying gridded emissions for ne0CONUSne30x8 simulations. Primarily supports sensitivity experiments that perturb emission magnitudes or temporal distributions.

See also [`functions/NEIfileDatetime_shift_BYweekofday.py`](../../../functions/NEIfileDatetime_shift_BYweekofday.py) for a utility that remaps NEI emission file timestamps by day-of-week (used to apply 2017 NEI emissions to other simulation years).

---

## Files

| File | Description |
|------|-------------|
| `ne30CONUSne30x8_scale_emissions_by_region_c20260212.py` | Scale gridded emissions within a specified geographic region by a multiplicative factor. Works on 0.15°×0.15° conservative-regridded emission files. Template for any regional emission perturbation experiment. |
| `ne30CONUSne30x8_scale_emissions_by_Canadaregion_c20260212.py` | Same as above, applied specifically to Canadian source regions (uses `grid_files/mask_CanadaProvinces_ne0CONUS_ne30x8.nc`). Designed for 2023 Canadian wildfire biomass burning sensitivity experiments. |

---

## Typical Use

```python
# Scale Canadian BB emissions by 15% in a specified province region
# Input:  conservatively regridded 0.15° emission files
# Output: modified emission NetCDF files at same resolution
# Mask:   grid_files/mask_CanadaProvinces_ne0CONUS_ne30x8.nc
```

For full anthropogenic emissions preprocessing (merging NEI 2017 with CAMS), see the `ACP_MUSICANEI_scripts/anthroemis_MergeNEI2017_to_CAMS/` directory in the offline script archive.
