# MUSICAv0 Workflows

Modular workflows for regional atmospheric chemistry experiment design, diagnostics, and model–observation evaluation using **MUSICAv0** with the **ne0CONUSne30x8 variable-resolution grid** (~14 km refinement over the contiguous United States).

This repository organizes emissions processing, model postprocessing, and evaluation tools in a grid-specific, reproducible structure designed for regional air quality studies.

---

## Scientific Scope

These workflows support:

- Emission sensitivity experiment design (e.g., anthropogenic and biomass-burning perturbations)
- Regional masking and grid handling for unstructured `ncol` grids
- Postprocessing of model output (column integration, diurnal analysis, regional aggregation)
- Model–observation comparison (surface networks and satellite products)
- Visualization of spatial and temporal diagnostics

The structure reflects the full lifecycle of a modeling study:

**experiment setup → model output processing → evaluation → interpretation**

---

## Repository Structure
### `grid_information/`
Grid-specific configuration for ne0CONUSne30x8, including masks, SCRIP grid handling, region definitions, and helper utilities for unstructured grids.

### `emissions/`
Tools for emissions preprocessing and sensitivity experiment generation, including regional masking, scaling, and temporal restructuring.

### `postprocessing/`
Diagnostics and analysis tools for model output, including column calculations, diurnal compositing, anomaly detection, and regional statistics.

### `model_evaluation/`
Model–observation comparison utilities, including bias metrics, correlation analysis, and stratified evaluation (e.g., urban–rural contrasts).

### `visualization/`
Reusable plotting utilities for maps, time series, and vertical diagnostics.

### `utilities/`
Shared helper functions used across modules (e.g., local time conversion, pressure interpolation, statistical helpers).

---

## Grid Configuration

All workflows under `ne0CONUSne30x8/` assume:

- Variable-resolution CAM/MUSICA grid
- Unstructured `ncol` indexing
- ~14 km refinement over CONUS with coarser resolution elsewhere
- Region-based masking defined on the refined grid

Adaptation would be required for uniform-resolution or alternative refinement configurations.

---

## Dependencies

Typical Python environment:

- Python ≥ 3.9  
- xarray  
- numpy  
- pandas  
- matplotlib  
- scipy  
- geopandas (for mask generation)

An example `environment.yml` can be provided for reproducibility.

---

## Notes on Data

This repository does **not** contain large model output or emissions files.  
Scripts assume user-supplied input paths and are designed to operate on externally stored datasets.

---

## Author

Madankui Tao  
Postdoctoral Associate, MIT EAPS  
Atmospheric chemistry modeling and satellite data analysis
