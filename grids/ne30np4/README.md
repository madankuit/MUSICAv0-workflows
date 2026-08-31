# grids/ne30np4/

Scripts and grid files for the **ne30np4** global spectral-element grid — the standard-resolution CESM/MUSICA global configuration (~1° uniform resolution, 48,602 `ncol` columns).

This directory is structured in parallel with `../ne0CONUSne30x8/` and follows the same subdirectory layout.

---

## Grid description

| Property | Value |
|----------|-------|
| Grid name | `ne30np4` |
| Nominal resolution | ~1° global |
| Column count (`ncol`) | 48,602 |
| Typical use | Global MUSICA simulations, boundary conditions for ne0CONUS runs |

---

## Directory layout

```
ne30np4/
├── grid_files/          # SCRIP files, region masks, regridding weight files
├── regridding/          # Scripts to regrid ne30np4 output to lat-lon grids
├── postprocessing/      # VCD calculation, file merging, surface extraction
├── model_evaluation/    # AQS/SLAMS matching, statistical evaluation
├── satellite_comparison/# TROPOMI AK application and bias diagnostics
├── plotting/            # Map visualization scripts
└── emissions/           # Emission scaling / preprocessing scripts
```

---

## Grid files

The ne30np4 SCRIP grid and CONUS masks **ship with this repository** under
[`grid_files/`](grid_files/) — see its README for the full list. Reference them
through `config/paths.py`, never by hard-coded path:

| Config constant | File |
|-----------------|------|
| `SCRIP_NE30NP4` | `ne30np4_091226_pentagons.nc` |
| `MASK_NE30NP4_CONUS_80KM` | `ne30np4_091226_pentagons_CONUSlandMaskedFalse_80kmBuffer.nc` |
| `MASK_NE30NP4_LOWER48_50KM` | `ne30np4_091226_pentagons_Lower48StatesCoastal50kmMaskedFalse.nc` |

Large derived ESMF weight files stay on the cluster and resolve through the same
config (`WEIGHTS_NE30_TO_1X1`, `WEIGHTS_FV09X125_TO_NE30`).

---

## Workflow scripts

The ne30np4 workflow currently in use is the CONUS background-ozone project,
which lives at the repository root: [`CONUSBGO3/`](../../CONUSBGO3/).
