# plotting/

Map visualization scripts for ne0CONUSne30x8 simulation outputs. These are example scripts showing how to use `functions/Plot_2D.py` for CONUS-domain diagnostics.

For the core plotting class and colormap utilities, see [`functions/Plot_2D.py`](../../../functions/Plot_2D.py).

---

## Files

| File | Description |
|------|-------------|
| `Plot_CompMaps.py` | Generate side-by-side or difference maps comparing two model runs (e.g., base case vs. sensitivity experiment). Produces absolute and relative difference panels across CONUS. Exports to PDF. |

---

## Usage Pattern

```python
from Plot_2D import Plot_2D
import xarray as xr

ds = xr.open_dataset('my_model_output.nc')
scrip_file = '../grid_files/ne0CONUS_ne30x8_np4_SCRIP.nc'

p = Plot_2D(scrip_file=scrip_file)
p.plot_SE(ds['NO2'], title='Surface NO₂', cmap='YlOrRd', vmin=0, vmax=50)
```

---

## Adding New Plots

To plot a new diagnostic:
1. Load your ne0CONUSne30x8 output (shape: `ncol` or `time × ncol`)
2. Instantiate `Plot_2D` with the SCRIP file
3. Call `.plot_SE()` for unstructured SE data or `.plot_FV()` for regridded lat/lon data
4. Customize projection, colormap, and region as needed
