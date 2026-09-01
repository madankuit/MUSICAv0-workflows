# MUSICAv0 Workflows

Analysis and postprocessing code for **MUSICAv0** (CESM2 with CAM-chem in spectral-element mode), with a focus on regional air quality studies over the continental United States.

Scripts here cover the full modeling lifecycle:

**experiment setup → emissions processing → model output postprocessing → evaluation → visualization**

---

## Guides

**New to spectral-element output? [Working with the unstructured grid (`ne0CONUSne30x8`)](docs/working_with_ne0CONUSne30x8.md)**
— grid conventions, reading, plotting, extracting at a point, VCDs, and using 3-D
fields as satellite a priori, each step pointing at the script here that does it.

---

## Repository Structure

```
MUSICAv0-workflows/
├── config/
│   └── paths.py                 # ★ single source of truth for EVERY path
├── check_paths.py               # verify those paths exist on your machine
├── docs/
│   └── working_with_ne0CONUSne30x8.md   # guide to the unstructured grid
├── functions/                   # General-purpose, grid-agnostic utility functions
├── CONUSBGO3/                   # Project: CONUS background ozone (ne30np4)
│   ├── experiment_setup/        # CESM namelists + input preparation
│   ├── emissions_perturbation/  # build the noAnthro / noBB emission inputs
│   ├── postprocessing/          # monitor matching + MDA8 extraction
│   ├── regridding/              # ne30 → 1° conservative regrid + gridded MDA8
│   └── analysis/                # scenario differences and figures
└── grids/
    ├── ne0CONUSne30x8/          # Variable-resolution grid (~14 km over CONUS)
    │   ├── grid_files/          # SCRIP grid, land/region masks (.nc)
    │   ├── regridding/          # Regrid SE outputs to regular lat/lon
    │   ├── emissions/           # Emissions preprocessing and sensitivity experiments
    │   ├── postprocessing/      # Model output processing (VCD, surface extraction, merging)
    │   ├── model_evaluation/    # Model–observation comparison (AQS, SLAMS, MDA8 O3)
    │   ├── satellite_comparison/# TROPOMI column comparison with averaging kernels
    │   └── plotting/            # Map visualization scripts
    ├── ne30np4/                 # Standard-resolution global grid (~1°)
    │   └── grid_files/          # SCRIP grid + CONUS land masks (.nc)
    └── target_grids/
        └── grid_files/          # Regular lat/lon ESMF destination grids (.nc)
```

### `functions/`
Reusable utility functions that are **grid-agnostic** or work across multiple SE grids (ne30, ne0CONUSne30x8, etc.). Import these from any analysis script. See [functions/README.md](functions/README.md).

### `grids/ne0CONUSne30x8/`
**Grid-specific** workflow scripts for the `ne0CONUSne30x8` variable-resolution CAM-SE grid (174,098 unstructured `ncol` columns; ~14 km over CONUS, ~111 km globally). Each subfolder targets a distinct stage of analysis. See [grids/ne0CONUSne30x8/README.md](grids/ne0CONUSne30x8/README.md).

### `grids/ne30np4/`
Grid files for the `ne30np4` standard-resolution global CAM-SE grid (48,602 `ncol` columns; ~111 km globally): the SCRIP description and the CONUS land masks, all shipped in the repo. Used for global simulations and as a lateral boundary condition source for ne0CONUS runs. The workflow that runs *on* this grid is [`CONUSBGO3/`](CONUSBGO3/README.md). See [grids/ne30np4/README.md](grids/ne30np4/README.md).

### `CONUSBGO3/`
**A worked example of extracting surface ozone** from a MUSICAv0 `ne30np4` run — hourly surface O₃ and daily MDA8, at monitor points and on a regular 1° grid — built around a set of emission-zeroing sensitivity runs. Useful as a template for point extraction on an unstructured grid, local-time MDA8, and conservative SE → lat/lon regridding. See [CONUSBGO3/README.md](CONUSBGO3/README.md).

### `grids/target_grids/`
Regular lat/lon grid descriptions used as the **destination** of ESMF regridding (1°, 0.15°, 0.1° CAMS). Small files, shipped in the repo. See [grids/target_grids/README.md](grids/target_grids/README.md).

---

## Grid Overview

| Grid | Columns (`ncol`) | Resolution | Primary use |
|------|-----------------|------------|-------------|
| `ne0CONUSne30x8` | 174,098 | ~14 km CONUS / ~111 km global | Regional air quality |
| `ne30np4` | 48,602 | ~111 km global | Global / LBC source runs |
| `f09` | — | 0.9°×1.25° FV | Finite-volume reference |

SCRIP and mask files for ne0CONUSne30x8 live in [grids/ne0CONUSne30x8/grid_files/](grids/ne0CONUSne30x8/grid_files/); ne30np4 grid files in [grids/ne30np4/grid_files/](grids/ne30np4/grid_files/).

---

## Getting Started

### 1. Clone and set up the environment

```bash
git clone https://github.com/madankuit/MUSICAv0-workflows.git
cd MUSICAv0-workflows
conda env create -f environment.yml
conda activate musica-workflows
```

### 2. Configure your paths

**Every path used anywhere in this repository is declared in one file:**

```
config/paths.py
```

No script hard-codes a path — they all import from here. The file separates two kinds:

- **Repo-relative** — SCRIP grids, land/region masks, ESMF destination grids and the
  shared `functions/` all ship with the repository and are addressed relative to the
  repo root. These work on any machine, in any checkout location, **with no
  configuration at all**.
- **Cluster data** — model archive, satellite data, emissions and processed output live
  outside the repo. These are built from the *runtime* user name, so the defaults are
  already correct if your data sits under `/net/fs09/d0/$USER`.

**Where new cluster output belongs.** Cluster paths split into two kinds, and the
difference matters when you add one:

| Kind | Goes under | Examples |
|------|-----------|----------|
| **This project's own data** | `MUSICA_PROJECT_ROOT` — `$DATA_ROOT/MUSICAv0-workflows/{data,figures}/` | `MUSICA_COLIDX_DIR` (AQS monitor → column-index CSVs) |
| **Shared, cross-project data** | the data-type folders — `CESM22/`, `Datasets/`, `ncar_copies/` | model archive, TROPOMI L2/L3 products, emission inventories |

Put **new** MUSICA-specific outputs under `MUSICA_PROJECT_ROOT`. Satellite and
emission inputs stay in the shared folders — other projects read them too, and
splitting related products apart (the NO₂ AK and AMF files, say) is how pairings
get broken.

A few existing outputs predate this convention and deliberately stay put, because
moving them would break paths for no scientific gain: `BGO3_ROOT`
(under `ProcessedData/`), `TROP_VCD_ROOT`/`TOTAL_VCD_ROOT`
and `REGRIDDED_2018_ROOT` (under `CESM22/`). Follow the rule above for anything
new rather than matching them.

If your data lives somewhere else, override it with an environment variable rather than
editing the file:

```bash
export MUSICA_ENV_DATA_ROOT=/net/fs09/d0/<someone>   # large data volume
export MUSICA_ENV_HOME_ROOT=/home/<someone>          # home / scripts volume
```

The other `MUSICA_ENV_*` overrides (individual leaf paths, and the provenance strings
written into output NetCDFs) are documented at the top of `config/paths.py`.

To check what actually exists on your machine:

```bash
python check_paths.py            # repo-shipped files and cluster data
python check_paths.py --repo     # only the files that ship with the repo
python check_paths.py --suggest  # also hunt for moved files/directories
```

Repo-shipped files should always be present; a miss there means an incomplete checkout.
Cluster misses are expected on a machine that does not hold the data.

### 3. Importing the config from a script

Scripts locate the config by walking up from their own location, so this works from any
subdirectory and contains no absolute path:

```python
import sys, pathlib
_ROOT = next(p for p in pathlib.Path(__file__).resolve().parents
             if (p / 'config' / 'paths.py').exists())
sys.path.insert(0, str(_ROOT))
import config                       # also puts functions/ on sys.path
from config.paths import ARCHIVE, SCRIP_NE30NP4
```

Because `import config` adds `functions/` to `sys.path`, shared utilities import by plain
module name:

```python
from SE_analysis import get_site_index
```

In notebooks the same bootstrap is used, walking up from `Path.cwd()` instead of `__file__`.

### 4. Typical analysis workflow (ne0CONUSne30x8)

```
1. Regrid SE output to lat/lon         →  grids/ne0CONUSne30x8/regridding/
2. Compute VCDs / surface extractions  →  grids/ne0CONUSne30x8/postprocessing/
3. Apply TROPOMI averaging kernels     →  grids/ne0CONUSne30x8/satellite_comparison/
4. Match with AQS/SLAMS observations   →  grids/ne0CONUSne30x8/model_evaluation/
5. Plot maps and diagnostics           →  grids/ne0CONUSne30x8/plotting/
```

Each subfolder has its own README with script-level descriptions.

> **Note (TROPOMI NO₂ averaging kernel).** In step 3, the kernel stored in the TROPOMI L2
> product is the *total-column* kernel and is converted to a tropospheric one via
> `AK_trop = AK_total × (AMF_total / AMF_trop)` (TROPOMI ATBD S5P-KNMI-L2-0005-RP). This is
> implemented and applied for **NO₂ only** — HCHO and CO use their stored kernel directly.
>
> It needs `AMF_total`/`AMF_trop` on the AK's grid, and **those regridded fields do not exist
> yet, so NO₂ cannot currently be run.** It is a small gap rather than a data problem: both are
> ordinary L2 variables (`air_mass_factor_total`, `air_mass_factor_troposphere`) sitting in the
> same `/PRODUCT/` group as `averaging_kernel`, so the fix is to re-run the existing L2→L3
> regrid with those two added to its variable list. See
> [`functions/README.md`](functions/README.md) and
> [`grids/ne0CONUSne30x8/satellite_comparison/README.md`](grids/ne0CONUSne30x8/satellite_comparison/README.md).

---

## Dependencies

[environment.yml](environment.yml) is the single, authoritative list — it is kept
in step with what the code actually imports, so don't maintain a second copy here.

```bash
conda env create -f environment.yml
conda activate musica-workflows
python -c "import xarray, geopandas, timezonefinder, pytz, seaborn; print('OK')"
```

Beyond the usual scientific stack (xarray / numpy / pandas / matplotlib / cartopy /
scipy / geopandas), one group is easy to overlook:

- **`timezonefinder` + `pytz`** — required by every EPA MDA8 pipeline, which converts
  UTC to local time *before* computing the 8-hour maxima. Without them the ozone
  postprocessing cannot run at all.
**`esmpy` is only needed to _generate_ ESMF weights** (the `gen_*_weights.py`
scripts). Applying weights that already exist is a sparse matrix multiply and needs
nothing beyond scipy — which is why the production regrid runs in a plain
environment while weight generation needs a dedicated one (`rootxesmf` on Svante).

---

## Resources

### The model

MUSICA is NSF NCAR's framework for representing chemistry across atmosphere
models. **MUSICAv0** — the part this repository works with — is the CAM-chem
configuration on a regionally-refined spectral-element grid, ~1° globally
refining to ~14 km over CONUS, available in CESM2.2.

| | |
|---|---|
| [MUSICA project page](https://www2.acom.ucar.edu/sections/musica) | the umbrella project and its components |
| [MUSICAv0](https://www2.acom.ucar.edu/sections/musica-v0) | the configuration used here |
| [MUSICA wiki](https://wiki.ucar.edu/display/MUSICA/MUSICA+Home) | **the documentation for setting up and running it** |
| [MUSICA projects](https://www2.acom.ucar.edu/sections/musica-projects) | community activity |
| [CAM-chem](https://www2.acom.ucar.edu/gcm/cam-chem) · [wiki](https://wiki.ucar.edu/display/camchem/Home) | the chemistry model MUSICAv0 configures |
| [`ESCOMP/CAM`](https://github.com/ESCOMP/CAM) | the atmosphere model source |

### Learning and analysing the output

| | |
|---|---|
| [`NCAR/musica-tutorial`](https://github.com/NCAR/musica-tutorial) | NCAR's tutorial notebooks for working with MUSICAv0 output — the closest upstream analogue to this repository |
| [`NCAR/MUSICA-Tools`](https://github.com/NCAR/MUSICA-Tools) | NCAR's input (IPT) and output (TPO) processing code for CAM-chem |

> **Attribution.** Two scripts here descend directly from the Nov 2021 MUSICA
> tutorial — `grids/ne0CONUSne30x8/regridding/Regrid_SE_015FV.py` (from
> `Rewrite_output.ipynb`) and `.../ConservativeRegrid_TROPOMITime_MUSICAOutputs.py`
> (from `ConservativeRegrid_MUSICAOutputs.ipynb`) — as their docstrings record.

### Evaluating against observations

**MELODIES MONET** is the NCAR–NOAA diagnostic framework for comparing model
output against surface, aircraft and satellite observations. It covers much of the
same ground as `model_evaluation/` here, in a model-agnostic package, and is worth
reaching for before writing new comparison code.

| | |
|---|---|
| [`NCAR/MELODIES-MONET`](https://github.com/NCAR/MELODIES-MONET) | the tool |
| [documentation](https://melodies-monet.readthedocs.io/en/stable/) | user and developer guides |
| [MELODIES tutorial 2024](https://www2.acom.ucar.edu/events/melodies-tutorial-2024) | tutorial materials |
| [NOAA CSL page](https://csl.noaa.gov/groups/csl4/modeldata/melodies-monet/) | the NOAA side of the project |

### Elsewhere in the MUSICA family

Chemistry software developed alongside MUSICAv0. Nothing here depends on it, but
it shares the name and is the natural next stop if you move from a configured
CESM run toward the chemistry itself.

| | |
|---|---|
| [`NCAR/musica`](https://github.com/NCAR/musica) | *Multi-Scale Infrastructure for Chemistry and Aerosols* — the MICM solver and TUV-x photolysis, with Python/Fortran/C++/JS interfaces |
| [`NCAR/music-box`](https://github.com/NCAR/music-box) | **MusicBox**, NCAR's box/column model running MICM chemistry standalone — [docs](https://ncar.github.io/music-box/branch/main/) · [ACOM page](https://www2.acom.ucar.edu/modeling/musicbox) · [wiki](https://wiki.ucar.edu/spaces/MusicBox/overview) |
| [`NCAR/music-box-interactive-desktop`](https://github.com/NCAR/music-box-interactive-desktop) | MusicBox as a desktop application |

---

## Data

Model output, emissions files, and observational datasets are **not** included in this
repository — they are located through `config/paths.py`.

Grid support files **are** included: SCRIP grids, CONUS land/region masks and the regular
lat/lon ESMF destination grids live under `grids/*/grid_files/` and need no configuration.
Large derived ESMF **weight** files (3.6 MB–165 MB) are the exception — they stay on the
cluster and are regenerated by the `gen_*_weights.py` scripts when missing.

---
