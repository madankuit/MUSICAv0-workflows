# CONUSBGO3 — CONUS Background Ozone with MUSICAv0 (ne30)

Emission-zeroing sensitivity experiments with **MUSICAv0** to estimate the
**U.S. background contribution to surface ozone** at a network of monitors
provided by the **Dan Jaffe** group (Univ. of Washington Bothell). For each
monitor we produce simulated hourly surface O₃ and daily **MDA8 O₃** for a
**BASE** run and two perturbation runs in which anthropogenic (**noAnthro**)
or biomass-burning (**noBB**) emissions are removed over CONUS, for the
ozone seasons **April–October of 2022 and 2023**.

The difference `BASE − perturbation` isolates the O₃ **produced from** CONUS
anthropogenic / fire emissions; what remains in the perturbation run is the
**background** (imported + biogenic + non-CONUS) contribution.

---

## Model configuration

| Item | Value |
|------|-------|
| Model | MUSICAv0 = CESM2.2, CAM-chem (MOZART-TS1) in spectral-element (SE) mode |
| Grid | `ne30np4` — quasi-uniform global ~111 km (48,602 `ncol` columns) |
| Compset | `FCnudged` (`f.e22.FCnudged.ne30_ne30_mg17`) |
| Meteorology | Nudged to **MERRA-2** (`u,v,T` coef 0.25), Nov 2021 – Dec 2023 |
| Anthropogenic emis | CAMS-GLOB-ANT v6.2 (monthly, `ne30np4`) |
| Biomass burning | QFED 2.6 / FINN (daily, `ne30np4`) |
| Aircraft / volcano / other | CMIP6 + CAMS v5.1 (time-shifted to approximate 2022–2023) |
| Lower boundary (GHG/ODS) | `LBC_17500116-25001216_CMIP6_SSP585_0p5degLat` |
| Output stream | `h2` hourly; surface level (`lev = -1`) used for O₃ |

> ⚠️ `ne30np4` here is the **standard global grid**, *not* the variable-resolution
> `ne0CONUSne30x8` used elsewhere in this repo. CONUS emissions are perturbed on
> the global grid using a CONUS land mask + 80 km buffer.

---

## Experiment matrix

All cases are named `f.e22.FCnudged.ne30_ne30_mg17.BGO3.<scenario>`, archived under
`$svante_archive` (see [Paths](#paths)).

| Label | Case suffix | Period | Description |
|-------|-------------|--------|-------------|
| — | `GetInitialForY202204.01` | Nov 2021 → Apr 2022 | Spin-up to generate the initial condition for 2022-04-01 |
| **BASE2022** | `BASEY20220401TY20230401` | 2022-04-01 → 2023-04-01 | Baseline, all emissions on |
| **BASE2023** | `BASEY20230401TY20231101` | 2023-04-01 → 2023-11-01 | Baseline, all emissions on |
| **noAnthro2022** | `noANTHROemisCONUS80kmBufferY20220401TY20221101` | 2022 season | CONUS anthropogenic emissions = 0 |
| **noAnthro2023** | `noANTHROemisCONUS80kmBufferY20230401TY20231101` | 2023 season | CONUS anthropogenic emissions = 0 |
| **noBB2022** | `noBBemisCONUS80kmBufferY20220401TY20221101` | 2022 season | CONUS biomass-burning emissions = 0 |
| **noBB2023** | `noBBemisCONUS80kmBufferY20230401TY20231101` | 2023 season | CONUS biomass-burning emissions = 0 |
| _(aux)_ | `noBBemisGlobalY20220401TY20231231` | 2022-04 → 2023-12 | Global BB removed — extra sensitivity check |

"CONUS 80 km buffer" = the lower-48 states polygon (US Census `cb_2018_us_state_500k`)
dilated by an 80 km buffer; emissions are zeroed only inside this mask.

---

## Directory layout

```
CONUSBGO3/
├── experiment_setup/          # CESM namelists + input-file preparation
│   ├── namelists/             # user_nl_cam for each production case
│   │   ├── BGO3.GetInitialForY202204.01_user_nl_cam
│   │   ├── BGO3.BASEY20220401TY20230401_user_nl_cam
│   │   ├── BGO3.noANTHROemisCONUS80kmBuffer_user_nl_cam
│   │   ├── BGO3.noBBemisCONUS80kmBuffer_user_nl_cam
│   │   ├── BGO3.noBBemisGlobal_user_nl_cam
│   │   ├── LBC_user_nl_cam
│   │   └── drafts/            # scratch / test namelist iterations (provenance only)
│   ├── ModifyOtherFiles_nlFiles_forY2022T2023.ipynb   # time-shift CMIP6 aircraft/other emis to ~2022–23; rewrite nl file lists
│   └── CheckLBC.ipynb                                  # inspect the GHG/ODS lower-boundary file
│
├── emissions_perturbation/    # build the noAnthro / noBB emission inputs
│   ├── RemoveCONUSEmis_setup_v3.ipynb   # ★ main: build CONUS+80km land mask, zero ANT & BB over CONUS, tally CONUS emission totals
│   ├── EditANT6.2files.ipynb            # clean NaNs in CAMS-GLOB-ANT v6.2 ne30np4 files
│   └── archive/                         # superseded v1/v2 of the mask/removal notebook
│
├── postprocessing/            # monitor matching + MDA8 extraction  (main deliverable)
│   ├── GetMatched_ne30_DanJaffe_GivenMonitors_ColumnIndex.py   # 1. map monitors → ne30 column index
│   ├── Merge_h2files_hourlysurfO3.py                           # 2. merge hourly surface O3 per case
│   ├── Extract_givenmonitorO3_hourly_dailyMDA8_toNetCDF.py     # 3. extract hourly + compute MDA8 → NetCDF
│   └── Extract_MDA8O3_givenMonitors.ipynb                      # development notebook for step 3
│
├── regridding/                # ne30 → 1° mass-conservative CONUS grid + gridded MDA8
│   ├── gen_ne30_to_1x1_weights.py             # one-time: ESMF conservative weights (rootxesmf env)
│   ├── Regrid_ne30_surfO3_to_1x1_conserve.py  # ★ regrid 6 cases → hourly gridded + unified MDA8
│   ├── check_gridded_MDA8_deliverable.py      # QC: coverage/NaN report + single-year maps
│   ├── plot_gridded_MDA8_scenarios.py         # BASE + (BASE−noAnthro)/(BASE−noBB) figure
│   ├── plot_gridded_MDA8_3scenarios.py        # BASE/noAnthro/noBB absolute, shared scale
│   ├── archive/                               # single-day method trials (conservative vs linear)
│   └── README.md
│
└── analysis/                  # scenario differences & figures
    ├── ResultAnalysis_v20250707.ipynb   # monthly-mean surface-O3 differences (BASE − noAnthro, BASE − noBB), across years/regions
    └── FireAnalysis.ipynb               # fire (BB) emission maps & diagnostics
```

---

## Workflow

### 1. Experiment setup — `experiment_setup/`
1. **Spin-up** with the `GetInitialForY202204.01` namelist to create the 2022-04-01 initial condition.
2. **Prepare input emissions** — `ModifyOtherFiles_nlFiles_forY2022T2023.ipynb` time-shifts CMIP6
   aircraft/other-surface emission files (+10 yr) and copies CAMS v5.1 2021 → 2023 to approximate the
   2022–2023 period, then rewrites the `srf_emis_specifier` / `ext_frc_specifier` lists in the namelists.
3. `CheckLBC.ipynb` verifies the SSP585 lower-boundary (CH₄, etc.) file.
4. Run **BASE** for each season using the `BASE*` namelists.

### 2. Emissions perturbation — `emissions_perturbation/`
- `RemoveCONUSEmis_setup_v3.ipynb` (★) downloads the US Census state shapefile, builds a CONUS
  land mask with an 80 km buffer on `ne30np4`, then writes perturbed emission files with
  anthropogenic (CAMS-GLOB-ANT v6.2) **or** biomass-burning (QFED/FINN) emissions set to zero
  over CONUS land. It also tallies the removed CONUS emission totals per species.
  - Mask files produced (data, on Svante — not in this repo):
    `ne30np4_091226_pentagons_CONUSlandMaskedFalse_{50,80}kmBuffer.nc`, stored alongside the
    SCRIP grid in `~/Scripts/CESM_analysis/functions/`.
- `EditANT6.2files.ipynb` strips residual NaNs from the CAMS-GLOB-ANT v6.2 `ne30np4` files.
- Run **noAnthro** and **noBB** for each season with the matching namelists.

### 3. Postprocessing — `postprocessing/`  (run in order)
1. **`GetMatched_ne30_DanJaffe_GivenMonitors_ColumnIndex.py`**
   Reads the monitor list `MonitorInfo/Lee_Jaffe_GAM_stats.csv` and, using
   `get_site_index` (from `functions/SE_analysis.py`) on the `ne30np4` SCRIP grid,
   finds the nearest model column for each AQS monitor.
   → `MonitorInfo/MatchedMonitors_ne30_ColIdx.csv`
2. **`Merge_h2files_hourlysurfO3.py`**
   Concatenates the hourly `h2` history files for MM-DD `0401`–`1101` and keeps the surface
   layer O₃ for each case.
   → `h2_surfO3_merged/<case>.cam.h2.surflev.O3.<start>T<end>.nc`
3. **`Extract_givenmonitorO3_hourly_dailyMDA8_toNetCDF.py`**
   For every matched monitor, extracts the hourly surface O₃ column, converts to local time
   (fixed summertime UTC offset via `timezonefinder`), and computes daily **MDA8 O₃** following
   EPA convention (8-h rolling mean requiring ≥6 valid hours; a day is valid if ≥13 of the 17
   windows ending 07–23 local time are present). Writes two CF-compliant NetCDFs per case.
   → `ForGivenMonitors/LocalTimeMDA8O3.<label>.GivenMonitors.<start>T<end>.nc`  (ppb, local-time dates)
   → `ForGivenMonitors/UTChourlyO3.<label>.GivenMonitors.<start>T<end>.nc`      (hourly, UTC)

`Extract_MDA8O3_givenMonitors.ipynb` is the exploratory notebook the `.py` was distilled from.

### 4. Regridding to a 1° CONUS grid — `regridding/`
Spatial companion to step 3: the same surface O₃ / MDA8, but for the **entire CONUS** on a
regular **1°×1°** grid via **mass-conservative** ESMF remap (ne30np4 → 1° FV). See
[regridding/README.md](regridding/README.md). Products land in
`…/DanJaffeMUSICAPostprocessing/Regridded1deg/`:
- 6 per-case hourly gridded surface O₃ files (`hourly/`, UTC, ppb);
- one **unified** daily-MDA8 file `MDA8O3(scenario={BASE,noAnthro,noBB}, time, lat, lon)`
  spanning Apr–Oct 2022 + 2023 — the shareable product.

The gridded MDA8 uses the identical local-time-adjusted EPA method as step 3 and was
validated against the point deliverable (nearest cell vs monitor: r = 0.994, bias −0.11 ppb).
Before sharing, run `check_gridded_MDA8_deliverable.py` for the per-(scenario, year)
coverage and NaN report.

### 5. Analysis — `analysis/`
- `ResultAnalysis_v20250707.ipynb` — monthly-mean surface-O₃ differences between scenarios
  (`noAnthro − BASE`, `noBB − BASE`) and across years, with regional maps.
- `FireAnalysis.ipynb` — fire emission magnitude maps and diagnostics.
- Figures written to `Figures/CESM_analysis/BGO3/`.

---

## Paths

Every path below is defined once, in [`config/paths.py`](../config/paths.py) — no
script hard-codes a location. Paths come in two kinds:

**Ships with the repository** (repo-relative, valid on any machine, no setup):

| Purpose | Config constant | Location in repo |
|---------|-----------------|------------------|
| `ne30np4` SCRIP grid | `SCRIP_NE30NP4` | `grids/ne30np4/grid_files/ne30np4_091226_pentagons.nc` |
| CONUS land mask (80 km buffer) | `MASK_NE30NP4_CONUS_80KM` | `grids/ne30np4/grid_files/…_CONUSlandMaskedFalse_80kmBuffer.nc` |
| Lower-48 coastal mask (50 km) | `MASK_NE30NP4_LOWER48_50KM` | `grids/ne30np4/grid_files/…_Lower48StatesCoastal50kmMaskedFalse.nc` |
| 1°×1° destination grid | `FV_GRIDINFO_1X1` | `grids/target_grids/grid_files/FV1x1grid_info_c20241105.nc` |
| Shared functions | `FUNCTIONS_DIR` | `functions/` |

**Lives on the cluster** (built from the runtime `$USER`; override with the
`MUSICA_ENV_*` environment variables documented in `config/paths.py`):

| Purpose | Config constant | Default location |
|---------|-----------------|------------------|
| Model archive (case output) | `ARCHIVE` | `$DATA_ROOT/CESM22/archive/` |
| Project processed-data root | `BGO3_ROOT` | `$DATA_ROOT/ProcessedData/DanJaffeMUSICAPostprocessing/` |
| ↳ Monitor list (input) | `BGO3_MONITOR_LIST` | `…/MonitorInfo/Lee_Jaffe_GAM_stats.csv` |
| ↳ Matched column index | `BGO3_MONITOR_COLIDX` | `…/MonitorInfo/MatchedMonitors_ne30_ColIdx.csv` |
| ↳ Merged hourly surface O₃ | `BGO3_MERGED_SURFO3_DIR` | `…/h2_surfO3_merged/` |
| ↳ Per-monitor deliverables | `BGO3_GIVEN_MONITORS_DIR` | `…/ForGivenMonitors/` |
| ↳ Gridded 1° deliverables | `BGO3_REGRIDDED_1DEG_DIR` | `…/Regridded1deg/` (unified MDA8 + `hourly/`) |
| ne30→1° conservative weights | `WEIGHTS_NE30_TO_1X1` | `$DATA_ROOT/CESM22/grids/ESMFmap_ne30np4_TO_1x1_conserve_c20260708.nc` |
| Emission source files | `BB_EMIS_NE30NP4_DIR`, `CAMS_V62_NE30NP4_DIR` | `$DATA_ROOT/ncar_copies/acom/MUSICA/emissions/` |
| Figures | `BGO3_FIGURE_DIR` | `$DATA_ROOT/Figures/CESM_analysis/BGO3/` |

`$DATA_ROOT` defaults to `/net/fs09/d0/$USER`. Verify what exists on your machine with:

```bash
python check_paths.py            # both groups
python check_paths.py --suggest  # also hunt for moved files
```

The case names and the Apr–Oct window are config constants too
(`BGO3_CASES`, `BGO3_CASE_LABELS`, `BGO3_START_MMDD`, `BGO3_END_MMDD`), so the
experiment matrix is defined in exactly one place.

The reference `user_nl_cam` namelists under `experiment_setup/namelists/` carry
`<MUSICA_ENV_DATA_ROOT>` / `<MUSICA_ENV_HOME_ROOT>` placeholders where absolute
input paths appeared; substitute your own before using them in a CESM case.

Data files (`.nc`, `.csv`, archive output) are **not** version-controlled — they live on Svante.

---

## Dependencies

- **Shared functions** in the repo's top-level [`functions/`](../functions/):
  - `SE_analysis.py` → `get_site_index` (nearest SE column for a lat/lon via SCRIP)
  - `func_ModelEval_statistical_tests.py`, `func_MUSICA_DefineRegion.py` (analysis)
  - Scripts reach `functions/` by importing `config`, which puts the directory on
    `sys.path`; no path juggling is needed. The SCRIP grid and CONUS-mask NetCDFs
    now ship in the repo as well (see [Paths](#paths)), so nothing outside the
    repository is required beyond the model output itself.
- **Python env**: the repo [`environment.yml`](../environment.yml)
  (`musica-workflows`) now covers everything these scripts import, including
  `timezonefinder` and `pytz` for the local-time MDA8 conversion — no extra
  manual installs.
  On Svante the point + gridded postprocessing runs in the `base` env; only the
  one-time conservative-weight generation
  (`regridding/gen_ne30_to_1x1_weights.py`) needs `esmpy`, for which a dedicated
  env exists (`rootxesmf`, esmpy 8.7). Applying the weights afterwards is a
  sparse mat-mul and does not need esmpy.

---

## Provenance

Scripts were developed under a personal `Scripts/CESM_analysis/ne30_CONUS_BackgroundO3/`
directory on Svante (2025) and reorganized here. The science is unchanged; what changed is
that every path now resolves through [`config/paths.py`](../config/paths.py) instead of being
hard-coded, the grid and mask NetCDFs ship with the repository, and provenance strings written
into output NetCDFs (`processed_by`, `contact`, the `<tag>` in file names) come from config
rather than being baked into the scripts.
