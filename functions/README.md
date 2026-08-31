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
| `Regridding_ESMF_v1.py` | General ESMF regridding framework. Supports bilinear and mass-conservative (patch) methods for SE → regular lat/lon. Called by grid-specific regridding scripts. |

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
| `func_VerticalRegrid_TROPOMIAK_toMUSICAlevs.py` | Regrid TROPOMI L2 averaging kernels (TM5 pressure grid) to MUSICA hybrid pressure levels. Output on 0.1°×0.1° lat/lon. Supports NO₂ and HCHO. **For NO₂ the stored total-column AK is converted to a tropospheric AK** (see note below); requires the regridded AMF directories `AMF_NO2_total_diri` / `AMF_NO2_trop_diri`. |
| `func_VerticalRegrid_TROPOMIAK_toMUSICAlevs_015latlon.py` | Same as above but outputs on 0.15°×0.15° lat/lon grid (matches conservative regridding output). NO₂ total→tropospheric AK conversion uses the regridded AMF files configured in `L3_015AMF_diri_dic`. |
| `func_VerticalRegrid_TotalCO_TROPOMIAK_toMUSICAlevs_015latlon.py` | Total CO column variant of the above. CO is a total-column product — its AK is used directly, no AMF conversion. |

> **NO₂ total-column → tropospheric averaging kernel.** The AK stored in the TROPOMI
> L2 `PRODUCT` group (and hence in the regridded AK files read here) is the
> *total-column* averaging kernel. The tropospheric AK used for the tropospheric-VCD
> comparison is obtained by scaling with the per-pixel air-mass-factor ratio:
>
> ```
> AK_trop(l) = AK_total(l) × (AMF_total / AMF_trop)
> ```
>
> (TROPOMI ATBD S5P-KNMI-L2-0005-RP / Product User Manual Eq. 4). The ratio is
> vertically uniform, so it multiplies every MUSICA layer. This matches the L2
> match/recalc convention (`tempo-v04-satellite-intercomparison`,
> `match_L2Scan._apply_ak_trop_conversion`). It is applied for **NO₂ only**; HCHO is a
> (tropospheric) total-column retrieval and uses its stored AK directly, as does CO.
> The conversion requires regridded `AMF_total` and `AMF_trop` fields on the same grid
> as the AK, produced by the upstream L2→L3 regrid step.

> **⚠️ Status (checked 2026-08-31): the regridded AMF inputs do not exist yet, so
> the NO₂ conversion cannot run.** The code is correct and in place; it has never
> had inputs.
>
> **This is a small gap, not a data problem.** `AMF_total` and `AMF_trop` are
> ordinary TROPOMI L2 variables — `air_mass_factor_total` and
> `air_mass_factor_troposphere` — and they sit in the **same `/PRODUCT/` group as
> `averaging_kernel`**, with the same per-pixel dimensions
> `(time, scanline, ground_pixel)`. Nothing extra has to be obtained or derived.
>
> What happened is simply that the L2→L3 regrid which produced the AK files
> carried the AK only: the regridded products hold `TROPOMI_NO2_AK` and no AMF, at
> both 0.15° and 0.1°. So **closing the gap means re-running that same regrid with
> the two AMF variables added to its variable list**, written on the AK's grid
> under the names this code expects (`TROPOMI_NO2_AMFtotal` /
> `TROPOMI_NO2_AMFtrop`).
>
> That is already demonstrated: the 0.05° CONUS `withAMF` product
> (`TROPOMI_005_WITHAMF_DIRS` in the config) carries them through as
> `total_air_mass_factor` / `tropospheric_air_mass_factor`, alongside the clear,
> cloudy and stratospheric AMFs and `tm5_tropopause_layer_index`. It is not a
> drop-in for this code — different grid, different variable names — but it shows
> the regrid change is a variable-list edit, not new work.
>
> **HCHO and CO are unaffected** — they use their stored AK directly.

### Utilities

| File | Description |
|------|-------------|
| `find_missing_files.py` | Identify missing hourly files within a datetime range. `find_missing_files_v1()` for NEI-format timestamps; `find_missing_files_v2()` for CAMS-format timestamps. |

---

## Notes

- `Regridding_ESMF_v1.py` requires `esmpy` (the Python interface to ESMF). Install via conda: `conda install -c conda-forge esmpy`.
- `Plot_2D.py` requires `cartopy`. Install via conda: `conda install -c conda-forge cartopy`.
- All functions accept `xarray.Dataset` or `numpy` arrays; output types match inputs where possible.
