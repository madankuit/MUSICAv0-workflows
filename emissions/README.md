# emissions/

Checks of the emission inputs read by MUSICAv0 runs. Paths come from
[`config/paths.py`](../config/paths.py); figures go to `EMISSIONS_CHECK_FIGURE_DIR`
(`$DATA_ROOT/Figures/EmissionsCheck/`).

| Script | What it does |
|---|---|
| `plot_bb_emissions_AprOct_BGO3.py` | QFED2.6 biomass-burning emissions (ne30np4) behind the CONUSBGO3 smoke-O3 attribution (`BASE − noBB`), Apr–Oct 2022 and 2023. For each species it maps the mean flux in the BASE input, the noBB input (CONUS land + 80 km buffer zeroed) and their difference, over North America. It also plots daily totals for the removed region vs fires left in both runs (Canada + Alaska, rest of the box) and writes the season totals to a CSV. |
| `run_plot_bb_emissions.sbatch` | SLURM wrapper (`base` env). From `emissions/`: `mkdir -p logs && sbatch run_plot_bb_emissions.sbatch [--species CO NO NO2 CH2O]` |

**Reading the smoke-O3 maps.** `BASE − noBB` only measures O₃ from fires *inside*
the zeroed region. O₃ from Canadian, Alaskan and Mexican fires is in both runs and
cancels, so this attribution misses transported smoke O₃ by construction.
