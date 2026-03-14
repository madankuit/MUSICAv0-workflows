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

SCRIP files used for ESMF regridding are stored under `grid_files/`.
See [`grid_files/README.md`](grid_files/README.md) for a list of available files.

Key files on Svante (paths in `svante_MUSICA_paths.py`):

| Variable | Filename |
|----------|----------|
| `SCRIP_ne30np4` | `ne30np4_091226_pentagons.nc` |
| `SCRIP_ne30np4_new` | `ne30np4_grid_c20241105.nc` |
| `Regridding_09x125weights_ne30` | `ESMFmap_0.9x1.25_ne30np4_cubit_conserve_cams.nc` |
