# Working with the MUSICAv0 unstructured grid (`ne0CONUSne30x8`)

A starting guide for collaborators who need to pull **hourly 3-D NO₂**, **surface
NO₂**, and **vertical profiles** out of a MUSICAv0 run on the variable-resolution
spectral-element grid, and use them as a priori for satellite retrievals.

It covers what the grid is, how to read and plot it, and which scripts in this
repository already do each step. Where something has *not* been done here, that
is said plainly rather than glossed over — see
[What is not covered](#what-is-not-covered-yet).

---

## 1. The run

```
case : f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.hourlyNEIemisCONUS.From20230101.01
path : $DATA_ROOT/CESM22/archive/<case>/atm/hist/
```

`$DATA_ROOT` is `/net/fs09/d0/<username>` on MIT Svante; in this repository it is
`ARCHIVE` in [`config/paths.py`](../config/paths.py), so no script hard-codes it.

**Coverage: 2023-01-01 → 2024-11-30**, which fully spans the requested
2023-12-01 → 2024-11-30 analysis year. 767 files in the archive.

### What already exists, and what you would generate

**Ready to use now — the raw model output.** Everything the request needs is in
the `h2` stream: hourly 3-D NO₂, surface NO₂, and the full vertical profile with
its meteorology. Nothing has to be run first to get at it.

**Not yet computed for this case.** There are no pre-computed VCDs, regridded
lat/lon fields, or monitor extractions for this run — those directories are empty
for it. The scripts to produce them are named in the sections below; running them
is expected to be the collaborator's own step, and the examples are there to work
from.

A useful reference point: the adjacent case
`…ne0CONUSne30x8_mt12.BASE.From20230801.01` **has** had tropospheric VCDs computed
with the very same script, so you can see the expected output format before
running anything:

```
$DATA_ROOT/CESM22/Calculated_MUSICA_VCD/TroposphericVCD/<that case>/
    <case>.cam.h2.TroposphericVCD.HCHO.20240501T20240531.nc
    -> variables: HCHO_TropVCD, lat, lon   dims: (time=744, ncol=174098)
```

One monthly file per species, still on the native `ncol` grid.

### History streams — use `h2`

| Stream | Frequency | Shape | Contents | Use it for |
|--------|-----------|-------|----------|------------|
| `h0` | monthly mean | 1 time/file, `(lev, ncol)` | ~900 variables | broad diagnostics |
| `h1` | daily | 31 times/file, `(ncol)` only | 2-D fields (`PM25_SRF`, `PRECT`, `T850`…) | quick surface/met checks |
| **`h2`** | **hourly** | **24 times/file, `(time, lev, ncol)`** | **3-D chemistry + full met** | **everything below** |

`h2` is the one you want: 700 files, one per day, 24 hourly steps each.

**3-D variables in `h2`** `(time, lev, ncol)`:

```
NO2  O3  CH2O  CO  PM25  CH3CN  HCN
T  Q  RELHUM  U  V  Z3  PMID  PDELDRY  MASS  M_dens  CLOUD
O3_Loss  O3_Prod  PAN_CHMP  r_NO2_OH  r_HO2_HO2  r_RO2_HO2
```

**2-D surface diagnostics** `(time, ncol)`: `NO2_SRF`, `NO_SRF`, `O3_SRF`,
`HNO3_SRF`, `PAN_SRF`, `OH_SRF`, `HO2_SRF` and others, plus `PS`, `TROP_P`.

> **This run carries its own meteorology.** `PS`, `T`, `PMID`, `TROP_P`, `Z3`,
> `PDELDRY` and the hybrid coefficients `hyam/hybm/hyai/hybi` are all in the `h2`
> file. That matters for picking the right VCD script — see
> [§6](#6-tropospheric-and-total-column-vcds). Some older cases did **not** save
> `PS` and needed an external met file; this one does not.

---

## 2. The grid, in one minute

`ne0CONUSne30x8` is a **spectral-element grid with regional refinement**. The
practical consequences:

- **There is no lat/lon axis.** Fields are `(time, lev, ncol)` — a single flat
  list of 174,098 columns. `lat` and `lon` are *variables* of length `ncol`, not
  dimensions. You cannot index by latitude, slice a box, or call `pcolormesh`
  directly.
- **Resolution varies across the domain.** Median cell width is ~17 km globally
  but ~11 km inside CONUS (min 5 km, max 139 km at the coarse global edge).
  76,237 of the 174,098 columns fall inside the CONUS box. Any area-weighted
  statistic must use the `area` variable — an unweighted mean over `ncol` is
  wrong, and silently so.
- **32 vertical levels**, `ilev` = 33 interfaces.

The grid geometry lives in a **SCRIP** file, which ships with this repository:

```python
from config.paths import SCRIP_NE0CONUSNE30X8   # grids/ne0CONUSne30x8/grid_files/
```

### Conventions that trip people up

| | Value | Note |
|---|---|---|
| `lon` | **0 – 360** | Not −180…180. Convert with `lon = (lon + 180) % 360 - 180`, or pass `360 + lon` when looking up a site with a negative longitude. |
| `lev` ordering | **TOA first** | `lev[0]` ≈ 3.6 hPa, `lev[-1]` ≈ 992.6 hPa. **Surface is index `-1`.** |
| `PS`, `PMID` | **Pa** | Not hPa. |
| `NO2`, `O3`, `CO`… | **mol/mol** | Multiply by 1e9 for ppbv. |

---

## 3. Reading the output

```python
import xarray as xr
from config.paths import ARCHIVE, case_hist_dir

case = 'f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.hourlyNEIemisCONUS.From20230101.01'
ds   = xr.open_dataset(case_hist_dir(case) / f'{case}.cam.h2.2024-06-15-00000.nc')

no2_3d   = ds['NO2']                    # (time=24, lev=32, ncol=174098), mol/mol
no2_surf = ds['NO2'].isel(lev=-1)       # (time, ncol) — surface, NOT lev=0
```

For multi-day work, open a set of files at once — but concatenate on `time`
explicitly, because `ncol` is shared and xarray will otherwise try to merge it:

```python
import glob
files = sorted(glob.glob(str(case_hist_dir(case) / f'{case}.cam.h2.2024-06-*.nc')))
ds = xr.open_mfdataset(files, combine='nested', concat_dim='time',
                       coords='minimal', compat='override')
```

Helper functions live in
[`functions/func_ReadMUSICAOutput.py`](../functions/func_ReadMUSICAOutput.py) —
`BasicRead_ds()` for a simple load, `Read_histh2_Runds()` for time/level-sliced
extraction.

> **A note on volume.** One `h2` file is 24 hours × 32 levels × 174,098 columns
> per 3-D variable. A full year is 365 of these. Subset variables and levels
> *before* concatenating, or you will exhaust memory. Starting with one week over
> one region (as suggested) is the right instinct.

---

## 4. Plotting an unstructured field

`pcolormesh` will not work — there is no 2-D mesh. You need the cell corners from
the SCRIP file and a polygon collection. That is what
[`functions/Plot_2D.py`](../functions/Plot_2D.py) does; it handles both
finite-volume and SE (including regionally-refined) output, with projections,
custom colormaps, log/symlog scaling and difference maps.

A worked example using it:
[`grids/ne0CONUSne30x8/plotting/Plot_CompMaps.py`](../grids/ne0CONUSne30x8/plotting/Plot_CompMaps.py).

**Quick-look alternative.** For a fast look without the polygon machinery, a
scatter of the column centres is usually enough and is honest about the grid:

```python
import matplotlib.pyplot as plt, cartopy.crs as ccrs
lon = (ds.lon.values + 180) % 360 - 180          # 0-360 -> -180..180
fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
ax.scatter(lon, ds.lat.values, c=ds['NO2'].isel(time=18, lev=-1) * 1e9,
           s=1, vmin=0, vmax=10, transform=ccrs.PlateCarree())
ax.coastlines(); ax.set_extent([-125, -66, 23, 50])
```

Use `Plot_2D.py` for anything publication-facing — the scatter misrepresents cell
area, which is exactly the thing that varies on this grid.

To cut a regional subset of columns, use `get_scrip_slice()` in
[`functions/SE_analysis.py`](../functions/SE_analysis.py).

---

## 5. Extracting at a point (surface NO₂ for column:surface relationships)

The nearest model column to a monitor is found with `get_site_index()` from
[`functions/SE_analysis.py`](../functions/SE_analysis.py). **Note the argument
order — longitude first**, and pass longitude in the model's 0–360 convention:

```python
from SE_analysis import get_site_index
from config.paths import SCRIP_NE0CONUSNE30X8

idx = get_site_index(site_lon=360 + monitor_lon,   # monitor_lon negative in CONUS
                     site_lat=monitor_lat,
                     scrip_file=str(SCRIP_NE0CONUSNE30X8))

no2_ppbv = ds['NO2'].isel(lev=-1, ncol=idx) * 1e9   # hourly surface NO2 at that site
```

It accounts for grid-cell boundaries rather than taking a naive nearest-centre,
which matters where resolution changes.

**Worked examples in this repository:**

| Script | What it does |
|--------|--------------|
| [`grids/ne0CONUSne30x8/model_evaluation/GetMatched_ne0CONUSne30x8_hourlyAQS_ColumnIndex.py`](../grids/ne0CONUSne30x8/model_evaluation/GetMatched_ne0CONUSne30x8_hourlyAQS_ColumnIndex.py) | Maps every AQS monitor to its `ncol` index and writes the lookup CSVs. **Run this once first** — everything else reuses the CSV. |
| [`grids/ne0CONUSne30x8/model_evaluation/MatchHourly_SLAMS_MUSICA_July2018_toCSV_v1.py`](../grids/ne0CONUSne30x8/model_evaluation/MatchHourly_SLAMS_MUSICA_July2018_toCSV_v1.py) | Hourly model-vs-observation matching at those monitors → CSV. |
| [`grids/ne0CONUSne30x8/model_evaluation/MatchDaily_SLAMS_MUSICA_July2018_toCSV.py`](../grids/ne0CONUSne30x8/model_evaluation/MatchDaily_SLAMS_MUSICA_July2018_toCSV.py) | Daily version of the same. |
| [`grids/ne0CONUSne30x8/postprocessing/Extract_Model_surfacelev.py`](../grids/ne0CONUSne30x8/postprocessing/Extract_Model_surfacelev.py) | Pulls the surface level out of `h2` for the whole domain. |
| [`grids/ne0CONUSne30x8/postprocessing/Merge_h2files.py`](../grids/ne0CONUSne30x8/postprocessing/Merge_h2files.py) | Concatenates many `h2` files into one time series. |

A simpler end-to-end pair, on the global `ne30np4` grid but structurally
identical and easier to read first:
[`CONUSBGO3/postprocessing/GetMatched_ne30_DanJaffe_GivenMonitors_ColumnIndex.py`](../CONUSBGO3/postprocessing/GetMatched_ne30_DanJaffe_GivenMonitors_ColumnIndex.py)
then
[`CONUSBGO3/postprocessing/Extract_givenmonitorO3_hourly_dailyMDA8_toNetCDF.py`](../CONUSBGO3/postprocessing/Extract_givenmonitorO3_hourly_dailyMDA8_toNetCDF.py)
(the latter also shows the UTC → local-time conversion, which must happen
*before* any daily statistic is computed).

> **`NO2` at `lev=-1` vs `NO2_SRF`.** `NO2.isel(lev=-1)` is the lowest model layer
> mid-point (~993 hPa, tens of metres thick); `NO2_SRF` is the model's own surface
> diagnostic. They are close but not identical. Pick one and state which — for
> column:surface ratios, the lowest-layer value is the usual choice.

---

## 6. Tropospheric and total column VCDs

Column density is the pressure-weighted integral of the mixing ratio:

```
VCD = Σ_k  X_k · ΔP_k / (g · m_air)
```

with the tropopause from `TROP_P` for the tropospheric column.

| Script | Use when |
|--------|----------|
| [`CalcTropVCD_h2_MUSICA_ne30CONUSne30x8_metwithinFile_c20260202.py`](../grids/ne0CONUSne30x8/postprocessing/CalcTropVCD_h2_MUSICA_ne30CONUSne30x8_metwithinFile_c20260202.py) | **← use this one for this run.** Meteorology is read from the `h2` file itself. |
| [`CalcTropVCD_h2_MUSICA_ne30CONUSne30x8_v2.py`](../grids/ne0CONUSne30x8/postprocessing/CalcTropVCD_h2_MUSICA_ne30CONUSne30x8_v2.py) | Older cases whose `h2` stream lacked `PS`; reads an external met file. |
| [`CalcTotalVCD_h2_MUSICA_ne30CONUSne30x8_v1.py`](../grids/ne0CONUSne30x8/postprocessing/CalcTotalVCD_h2_MUSICA_ne30CONUSne30x8_v1.py) | Total (not tropospheric) column — used for CO. |

Both write on the native `ncol` grid, so no regridding is needed if you are
staying unstructured.

### Getting to a regular lat/lon grid

If you need lat/lon (most satellite comparisons do), regrid **conservatively** —
linear interpolation between column centres does not conserve mass, which matters
for a column being compared against a retrieval:

- [`grids/ne0CONUSne30x8/regridding/ConservativeRegrid_TROPOMITime_MUSICAOutputs.py`](../grids/ne0CONUSne30x8/regridding/ConservativeRegrid_TROPOMITime_MUSICAOutputs.py)
  — SE → 0.15° lat/lon via offline ESMF weights.
- [`grids/ne0CONUSne30x8/regridding/Regrid_SE_015FV.py`](../grids/ne0CONUSne30x8/regridding/Regrid_SE_015FV.py)
  — the underlying regrid.
- Weight generation is illustrated for the global grid in
  [`CONUSBGO3/regridding/gen_ne30_to_1x1_weights.py`](../CONUSBGO3/regridding/gen_ne30_to_1x1_weights.py);
  a first-order conservative map is **not reversible**, so generate the direction
  you need.

---

## 7. 3-D NO₂ as satellite a priori

The model layer pressures you need are already in the file — `PMID` (Pa,
mid-layer) and `PS`, plus `hyai`/`hybi` if you want interfaces:

```
P_interface(k) = hyai(k)·P0 + hybi(k)·PS
```

(`P0` is not stored in `h2`; use the CAM standard 100000 Pa.)

So a per-column a priori profile is `NO2(time, :, ncol)` with `PMID(time, :, ncol)`
— no reconstruction needed.

### TROPOMI — examples exist here

The full chain is implemented:

1. [`functions/func_VerticalRegrid_TROPOMIAK_toMUSICAlevs_015latlon.py`](../functions/func_VerticalRegrid_TROPOMIAK_toMUSICAlevs_015latlon.py)
   — maps the TROPOMI averaging kernel from its TM5 pressure grid onto MUSICA
   levels (0.15° output; there is a 0.1° sibling).
2. [`grids/ne0CONUSne30x8/satellite_comparison/ConservativeRegrid_CalcTropVCD_approx1330LT_h2_MUSICA015deglatlon.py`](../grids/ne0CONUSne30x8/satellite_comparison/ConservativeRegrid_CalcTropVCD_approx1330LT_h2_MUSICA015deglatlon.py)
   — applies the AK and computes the AK-adjusted model column at ~13:30 LT.
3. [`grids/ne0CONUSne30x8/satellite_comparison/MaskTROPOMInan_in_MUSICA_VCD.py`](../grids/ne0CONUSne30x8/satellite_comparison/MaskTROPOMInan_in_MUSICA_VCD.py)
   — masks the model where the retrieval has no valid data, so the comparison is
   like-for-like.

**Two things to know before reusing this for NO₂:**

- The AK stored in the TROPOMI L2 `PRODUCT` group is the **total-column** kernel.
  For a tropospheric comparison it must be rescaled per pixel:
  `AK_trop = AK_total × (AMF_total / AMF_trop)` (TROPOMI ATBD
  S5P-KNMI-L2-0005-RP). This applies to **NO₂ only** — HCHO and CO use their
  stored kernel directly.
- That rescaling needs regridded `AMF_total`/`AMF_trop` fields, and **those do
  not currently exist** at the required grid. See the status note in
  [`functions/README.md`](../functions/README.md) for what is missing and what
  0.05° CONUS AMF data *is* available. HCHO and CO are unaffected.

---

## What is not covered yet

**TEMPO recalculation has not been done in this repository.** There is no TEMPO
code here — the satellite examples above are TROPOMI only. TEMPO differs in ways
that matter for NO₂, notably the **temperature correction factor**, so the
TROPOMI scripts should not be assumed to transfer unchanged.

For the exact TEMPO recalculation code, **Joe is the better reference** — he has
done that specific work. The TROPOMI examples here are a reasonable model for the
*structure* of the calculation (vertical regridding of the kernel, applying it to
the model profile, masking to valid retrieval pixels), but the TEMPO-specific
coefficients should come from him.

A fuller write-up of the satellite recalculation side will follow later.

---

## Suggested starting point

1. Open a single `h2` file, plot surface NO₂, confirm the 0–360 longitude and
   TOA-first level conventions look right in your plot.
2. Run the monitor-matching script once to build the `ncol` lookup CSV.
3. Pull one week (or one month) of hourly surface NO₂ over one region and get it
   to the rest of the team — enough to start the column:surface work without
   waiting on the full year.
4. Then scale to 2023-12-01 → 2024-11-30, subsetting variables and levels before
   concatenating.

## Environment

```bash
conda env create -f environment.yml
conda activate musica-workflows
```

See the repository [README](../README.md) for how paths are configured — every
path resolves through [`config/paths.py`](../config/paths.py), and
`python check_paths.py` reports what exists on your machine.
