# model_evaluation/

Compare ne0CONUSne30x8 model output against surface observations from the EPA AQS/SLAMS network. Covers grid-column matching, hourly/daily time series extraction, MDA8 O₃ computation, and export to CSV/NetCDF for statistical analysis.

For the statistical test functions themselves (Pearson/Spearman correlation, Wilcoxon test, RMA regression), see [`functions/func_ModelEval_statistical_tests.py`](../../../functions/func_ModelEval_statistical_tests.py).

---

## Files

### Grid-Column Matching

| File | Description |
|------|-------------|
| `func_MUSICAvsSLAMS_surf.py` | Core functions for ne0CONUSne30x8 model–observation comparison. `write_MonitorIDtoMUSICAcolidx_csvfile()` maps AQS monitor lat/lon to the nearest ne0CONUSne30x8 grid column using the SCRIP file. `plot_scatter_with_regression()` and `plot_difference_on_CONUSmap()` provide visualization. |
| `GetMatched_ne0CONUSne30x8_hourlyAQS_ColumnIndex.py` | Script to generate and save a CSV mapping AQS monitor IDs to ne0CONUSne30x8 `ncol` indices. Run once per AQS monitor set; output CSV is reused by matching scripts. |

### Time Series Matching and Export

| File | Description |
|------|-------------|
| `MatchHourly_SLAMS_MUSICA_July2018_toCSV_v1.py` | Extract model values at matched AQS monitor columns (hourly), align with SLAMS observations for July 2018, convert units, and export to CSV. |
| `MatchDaily_SLAMS_MUSICA_July2018_toCSV.py` | Same as above but for daily-mean values (h1 output). |

### O₃-Specific Extraction

| File | Description |
|------|-------------|
| `Extract_givenmonitorO3_hourly_dailyMDA8_toNetCDF.py` | Extract hourly surface O₃ at specified monitor locations from merged h2 files. Compute **daily MDA8 O₃** following EPA convention (max 8-hour rolling average, evaluated in local time). Output to NetCDF. Supports base and sensitivity scenarios. |

---

## Typical Workflow

```
1. GetMatched_ne0CONUSne30x8_hourlyAQS_ColumnIndex.py
   → monitors_to_ncol_index.csv

2. MatchHourly_SLAMS_MUSICA_July2018_toCSV_v1.py
   → matched_hourly_model_obs.csv

3. Use functions/func_ModelEval_statistical_tests.py
   → Pearson/Spearman correlation, bias, RMSE

4. For O3 background analysis:
   postprocessing/Merge_h2files_hourlysurfO3.py
   → Extract_givenmonitorO3_hourly_dailyMDA8_toNetCDF.py
   → daily MDA8 O3 at monitor sites (.nc)
```

---

## Notes

- `func_MUSICAvsSLAMS_surf.py` has hardcoded references to the ne0CONUSne30x8 SCRIP file path. Update `scrip_path` at the top of the file for your local data directory.
- Monitor-to-column matching is done using great-circle distance on the unstructured SE grid; column boundaries from the SCRIP file are used to verify containment.
- MDA8 computation follows 40 CFR Part 50 (8-hour rolling max, requiring ≥ 18 valid hourly values per day), applied in local standard time.
