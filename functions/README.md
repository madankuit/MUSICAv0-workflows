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

> **⚠️ Status (checked 2026-08-31): the AMF inputs do not exist yet, so the NO₂
> conversion cannot run.** The code is correct and in place; it has never had
> inputs. On the reference cluster:
> - 0.15°: the AK directory holds only `..._AK_Global_Regrid015deg_*.nc` files —
>   **no** `AMFtotal`, **no** `AMFtrop`;
> - 0.1°: the `TROPOMI_NO2AMFtotal_01deg` / `TROPOMI_NO2AMFtrop_01deg` directories
>   are absent (`check_paths.py` reports these two as MISSING).
>
> Attempting NO₂ therefore fails: the 0.1° entry point raises `ValueError` if the
> AMF directories are not supplied, and the 0.15° one raises `FileNotFoundError`.
> **HCHO and CO are unaffected** — they use their stored AK directly.
>
> **AMF fields that DO exist.** A 0.05° CONUS product with air-mass factors has
> already been produced — 34 daily files per species — and is addressed in the
> config as `TROPOMI_005_WITHAMF_DIRS`:
>
> ```
> $DATA_ROOT/Datasets/S5P_L2__NO2____HiR_2/withAMF_CONUS_005/
> $DATA_ROOT/Datasets/S5P_L2__HCHO___HiR_2/withAMF_CONUS_005/
> ```
>
> It is **not a drop-in** for the loaders here, which is why NO₂ still cannot run
> unmodified. The differences:
>
> | | the 0.05° product | what the loaders expect |
> |---|---|---|
> | grid | 0.05° CONUS | 0.15° (or 0.1°) global |
> | layout | one file per day, all AMFs (trop/strat/total/clear/cloud) | separate `AMFtotal_…` / `AMFtrop_…` files |
> | variables | TROPOMI L2 AMF names | `TROPOMI_NO2_AMFtotal` / `TROPOMI_NO2_AMFtrop` |
>
> So closing the gap needs either a 0.15°-global AMF regrid emitting those names,
> or a loader adapted to read the 0.05° CONUS files. **Neither has been done** —
> this note exists so the available data can be found, not because a path is
> wired up. The producing script is
> `TROPOMI_RegridL2HiR_UseOPeNDAP_MultiYear_005LatLon_{NO2,HCHO}_withAMF.py`,
> which lives outside this repository.


### Utilities

| File | Description |
|------|-------------|
| `find_missing_files.py` | Identify missing hourly files within a datetime range. `find_missing_files_v1()` for NEI-format timestamps; `find_missing_files_v2()` for CAMS-format timestamps. |

---

## Notes

- `Regridding_ESMF_v1.py` requires `esmpy` (the Python interface to ESMF). Install via conda: `conda install -c conda-forge esmpy`.
- `Plot_2D.py` requires `cartopy`. Install via conda: `conda install -c conda-forge cartopy`.
- All functions accept `xarray.Dataset` or `numpy` arrays; output types match inputs where possible.
