# satellite_comparison/

Compare ne0CONUSne30x8 model columns against TROPOMI L2 satellite retrievals. Applies TROPOMI averaging kernels to model profiles, computes bias and correlation metrics, and masks invalid satellite pixels consistently across datasets.

All scripts operate on **0.15°×0.15° conservatively regridded** model output (produced in `regridding/`).

---

## Files

### Averaging Kernel Application (VCD at TROPOMI Overpass Time)

| File | Species | Description |
|------|---------|-------------|
| `ConservativeRegrid_CalcTropVCD_approx1330LT_h2_MUSICA015deglatlon.py` | HCHO, NO₂ | Apply TROPOMI tropospheric averaging kernels to regridded model profiles at ~1:30 PM local time. Computes kernel-weighted tropospheric VCD on the 0.15° grid. **NO₂:** the AK is converted from total-column to tropospheric (`AK_trop = AK_total × AMF_total / AMF_trop`) inside `func_VerticalRegrid_TROPOMIAK_toMUSICAlevs_015latlon.py`; HCHO uses its stored AK directly. |
| `ConservativeRegrid_CalcTotalVCDCO_approx1330LT_h2_MUSICA015deglatlon.py` | CO | Same approach for total CO column using TROPOMI CO averaging kernels. |

### Data Preprocessing

| File | Description |
|------|-------------|
| `MaskTROPOMInan_in_MUSICA_VCD.py` | Mask model VCD pixels where corresponding TROPOMI pixels are NaN (cloud fraction, quality flag filters). Ensures consistent spatial coverage between model and satellite for fair comparison. |

### Metrics and Output

| File | Description |
|------|-------------|
| `WriteNCfile_ModelBiasToTROPOMI.TemporalCorrelation_CONUS.py` | Compute Spearman correlation, NMBE (normalized mean bias error), and NRMSE (normalized RMSE) between model and TROPOMI VCDs by US region and across CONUS. Outputs results to NetCDF. |

---

## Dependencies

- TROPOMI averaging kernels: from TROPOMI L2 offline products, mapped to 0.15°×0.15° (see `functions/func_VerticalRegrid_TROPOMIAK_toMUSICAlevs_015latlon.py`)
- TROPOMI AMF fields (NO₂ only): regridded `AMF_total` and `AMF_trop` on the same 0.15° grid as the AK, used to convert the stored total-column AK to a tropospheric AK. Configure via `L3_015AMF_diri_dic` in `func_VerticalRegrid_TROPOMIAK_toMUSICAlevs_015latlon.py`. Produced by the upstream L2→L3 regrid step alongside the AK.
  **⚠️ These regridded files do not exist yet — NO₂ cannot be run until they do.**
  This is a small gap: `air_mass_factor_total` and `air_mass_factor_troposphere`
  are standard TROPOMI L2 variables in the same `/PRODUCT/` group as
  `averaging_kernel`, so the fix is to re-run the existing L2→L3 regrid with those
  two added to its variable list. See
  [`functions/README.md`](../../../functions/README.md) for details.
  See the status note in [`functions/README.md`](../../../functions/README.md) for
  what is missing and why the existing `_withAMF` regrid is not a drop-in
  replacement. HCHO and CO are unaffected.
- Model VCDs: from `postprocessing/CalcTropVCD_h2_MUSICA_ne30CONUSne30x8_*.py`
- Regridded 0.15° model fields: from `regridding/ConservativeRegrid_*.py`

## Workflow

```
regridding/ConservativeRegrid_TROPOMITime_MUSICAOutputs.py
    ↓
ConservativeRegrid_CalcTropVCD_approx1330LT_h2_MUSICA015deglatlon.py
    ↓
MaskTROPOMInan_in_MUSICA_VCD.py
    ↓
WriteNCfile_ModelBiasToTROPOMI.TemporalCorrelation_CONUS.py
    ↓  regional bias / correlation metrics (.nc)
```

---

## Notes

- Averaging kernels are applied following the standard formulation: VCD_AK = Σ_k AK_k × x_model_k × Δp_k / g
- **NO₂ AK is tropospheric, not total-column.** The TROPOMI L2 `PRODUCT` AK is the total-column averaging kernel; for the tropospheric-VCD comparison it is scaled per pixel by `AK_trop = AK_total × (AMF_total / AMF_trop)` (TROPOMI ATBD S5P-KNMI-L2-0005-RP), matching the L2 match/recalc convention. HCHO (tropospheric total-column retrieval) and CO (total column) use their stored AK directly.
- The ~1:30 PM LT overpass time approximation uses UTC offsets from `functions/func_MUSICA_DefineRegion.py`
- TROPOMI quality flags: use `qa_value ≥ 0.75` (recommended by the TROPOMI L2 ATBD) for NO₂ and HCHO
