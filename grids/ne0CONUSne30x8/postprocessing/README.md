# postprocessing/

Process ne0CONUSne30x8 model history (h1/h2) files: compute column-integrated quantities, merge surface-layer outputs, and extract site-specific time series.

---

## Files

### Vertical Column Density (VCD) Calculations

These scripts compute tropospheric or total vertical column densities from hourly h2 history files. VCDs are integrated using hybrid pressure coordinates and model-internal meteorology.

| File | Species | Met source | Description |
|------|---------|-----------|-------------|
| `CalcTropVCD_h2_MUSICA_ne30CONUSne30x8_v2.py` | HCHO, NO₂ | External met file | Compute hourly tropospheric VCD using a separate meteorology file. Earlier approach (v2, July 2018 simulations). |
| `CalcTropVCD_h2_MUSICA_ne30CONUSne30x8_metwithinFile_c20260202.py` | HCHO, NO₂ | Within h2 file | Compute hourly tropospheric VCD reading met fields directly from h2 output. Updated approach for 2023–2024 simulations. **Preferred for new runs.** |
| `CalcTotalVCD_h2_MUSICA_ne30CONUSne30x8_v1.py` | CO, HCHO, NO₂ | External met file | Compute hourly **total** (full-column) VCD. |

### Surface-Layer Extraction and Merging

| File | Description |
|------|-------------|
| `Merge_h2files.py` | Merge multiple h2 hourly files across dates; extract surface layer (`lev=31`) for specified variables. Used for AQS/SLAMS model–observation matching. |
| `Merge_h2files_hourlysurfO3.py` | Merge h2 files for surface O₃ across multi-year date ranges. Supports base and sensitivity simulations (e.g., noBB, noAnthro). Designed for background O₃ analysis. |
| `Extract_Model_surfacelev.py` | Extract surface-level model variables from h2 files for a given set of dates and output to NetCDF. |

### Site Time Series

| File | Description |
|------|-------------|
| `CESM_WriteToHourlyDataframe.py` | Extract all CESM output variables at a specific location (site index) to a pandas DataFrame / CSV, in Eastern Daylight Time. |

---

## Notes

- VCD scripts output per-date NetCDF files in `ncol` coordinates (unstructured). Feed these into `satellite_comparison/` for TROPOMI averaging kernel application, or into `model_evaluation/` for surface comparisons.
- `Merge_h2files_hourlysurfO3.py` outputs are designed as input to `model_evaluation/Extract_givenmonitorO3_hourly_dailyMDA8_toNetCDF.py`.
- All scripts assume the `functions/func_ReadMUSICAOutput.py` reader is available on the Python path.
