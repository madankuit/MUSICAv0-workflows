# CONUSBGO3/regridding/

Mass-conservative regrid of MUSICAv0 **ne30np4 surface O3** to a regular **1°×1°
CONUS grid**, and derivation of **daily MDA8 O3**, for the BASE / noAnthro / noBB
experiments (2022 & 2023, Apr–Oct). Companion to [../postprocessing/](../postprocessing/)
(which extracts the same quantities at point monitors).

## Method

**Horizontal regrid** — ESMF **first-order conservative** remap (`RegridMethod.CONSERVE`),
ne30np4 → 1° FV, identical in spirit to the TROPOMI-overpass conservative regrid in
`grids/ne0CONUSne30x8/regridding/`. Applied as a sparse mat-mul of the offline weights.

**MDA8** — mirrors the point pipeline
[`Extract_givenmonitorO3_hourly_dailyMDA8_toNetCDF.py`](../postprocessing/Extract_givenmonitorO3_hourly_dailyMDA8_toNetCDF.py):
the hourly series is shifted from **UTC to local time first** (per-grid-cell summertime/DST
UTC offset via `timezonefinder`), **then** the EPA MDA8 is computed — 8-h rolling means
(≥6 valid hours), daily max over the 17 windows ending 07–23 local, day valid if ≥13/17.

## Grid files

| Role | Config constant | Where it lives |
|------|-----------------|----------------|
| Source SCRIP (SE) | `SCRIP_NE30NP4` | in the repo — `grids/ne30np4/grid_files/` |
| Destination 1° FV grid | `FV_GRIDINFO_1X1` | in the repo — `grids/target_grids/grid_files/` |
| **Conservative weights** ne30np4→1° | `WEIGHTS_NE30_TO_1X1` | on the cluster — `$DATA_ROOT/CESM22/grids/` (4 MB, derived) |

> The pre-existing `ESMFmap_0.9x1.25_ne30np4_cubit_conserve_cams.nc` from the emissions
> work is **FV→ne30** (the emissions direction) and cannot be reused for SE→FV; the
> ne30→1° weights above were generated fresh (conservative maps are not reversible).

## Scripts

| Script | Env | Purpose |
|--------|-----|---------|
| `gen_ne30_to_1x1_weights.py` | `rootxesmf` (esmpy 8.7) | One-time: generate + verify the ne30np4→1° conservative weight file. |
| `Regrid_ne30_surfO3_to_1x1_conserve.py` | `base` | Regrid the 6 merged surface-O3 files → hourly gridded + unified MDA8. Weight application is a sparse mat-mul (no esmpy needed). |
| `check_gridded_MDA8_deliverable.py` | `base` | **QC before sharing**: per-(scenario, year) day counts and NaN fractions, the Oct-31 local-time edge case, plus single-year BASE / BASE−noAnthro / BASE−noBB maps. |
| `plot_gridded_MDA8_scenarios.py` | `base` | 3-panel CONUS figure from the unified MDA8: (a) BASE, (b) BASE−noAnthro, (c) BASE−noBB, seasonal-mean over Apr–Oct 2022–2023. Independent diverging scales. |
| `plot_gridded_MDA8_3scenarios.py` | `base` | 3-panel CONUS figure of BASE / noAnthro / noBB as **absolute** seasonal means on one shared colour scale — the scenarios side by side rather than their differences. |

Run order: `gen_…_weights` (once) → `Regrid_…_conserve` → `check_…_deliverable`
→ either plotting script. All three consumers locate the unified MDA8 file by
glob (`bgo3_unified_mda8_glob()`), so none hard-codes its name or tag.

[`archive/`](archive/) holds the two single-day method-development trials — the
conservative one that was adopted and the linear-interpolation one that was
rejected. See its README for why.

## Outputs (Svante, `…/DanJaffeMUSICAPostprocessing/Regridded1deg/`)

- `hourly/CONUS1x1_UTChourlySurfO3.<label>.<start>T<end>.<tag>_c<YMD>.nc` — 6 files,
  hourly surface O3 `(time, lat, lon)` in **ppb**, **UTC** timestamps.
- `MUSICAv0_ne30_CONUS1x1_MDA8O3_BGO3_2022-2023_AprOct_<tag>_c<YMD>.nc` — **the shareable
  product**: `MDA8O3(scenario, time, lat, lon)` in ppb; `scenario ∈ {BASE, noAnthro, noBB}`,
  `time` = 428 daily local-time dates (Apr–Oct 2022 + Apr–Oct 2023), grid lat 24–50 °N,
  lon −125…−66 °W.

`<tag>` is `AUTHOR_TAG` from `config/paths.py` (defaults to the runtime `$USER`;
set `MUSICA_ENV_AUTHOR_TAG` to pin it). The plotting script locates the unified
file by glob, so it finds it whatever tag it carries.

Every file carries provenance in its global attributes (`processed_by`, `processing_date`,
`machine`, `source`, `horizontal_regrid`, `regrid_weight_file`, `mda8_method`, `local_time`,
`case_names`). Data stay on Svante — not version-controlled.

## Validation

Gridded MDA8 sampled at the nearest cell to each of the 781 monitors, vs the point
`LocalTimeMDA8O3` deliverable (BASE 2022): **mean bias −0.11 ppb, median |bias| 0.38 ppb,
median correlation 0.994** (98 % of sites r>0.9) — confirming the local-time-adjusted EPA
MDA8 is applied consistently with the point pipeline.

**Note:** ~24 of 1620 cells (Southern California; a few in ID/OR) show MDA8 > 150 ppb.
These reflect the **model's** surface-O3 extremes (present in the point files too — the
gridded/point agreement is r≈0.99), not a regrid artifact; the bulk field is realistic
(p99 ≈ 77 ppb). A single unrecoverable edge exists — Oct 31 in western time zones for
BASE-2023 and noBB-2022 is NaN (those merged inputs stop at Nov 1 00 UTC, so the last local
day has <13 valid 8-h windows and is correctly dropped) — 0.2 % of values.

## Figures

All figures land in `Figures/CESM_analysis/BGO3/regrid_trial/`.

`plot_gridded_MDA8_scenarios.py` → `MDA8_BASE_and_diffs_seasonmean_2022-2023.png`:
three CONUS panels — **(a) BASE**, **(b) BASE − noAnthro**, **(c) BASE − noBB** — as the
seasonal-mean MDA8 over Apr–Oct 2022–2023. Panels (b)/(c) use independent diverging scales
because the anthropogenic contribution (~30 ppb) is ~8× the biomass-burning contribution (~4 ppb).

`plot_gridded_MDA8_3scenarios.py` → `MDA8_3scenarios_seasonmean_2022-2023.png`:
the same period, but **BASE / noAnthro / noBB as absolute fields** on a single shared
colour scale (2nd–98th percentile of all three together), so the scenarios can be read
against one another directly instead of through their differences.

`check_gridded_MDA8_deliverable.py` → `DELIV_BASE_MDA8_seasonmean_<year>.png` and
`DELIV_BASE_minus_{noAnthro,noBB}_MDA8_seasonmean_<year>.png`: single-year versions used
as a sanity check on the deliverable, alongside the coverage report it prints.
