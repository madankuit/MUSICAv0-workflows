#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_bb_emissions_AprOct_BGO3.py

Biomass-burning (QFED2.6/FINN, ne30np4) emissions behind the CONUSBGO3 smoke-O3
attribution (BASE - noBB), over the same window as the smoke-O3 maps:
Apr 1 - Oct 31 of 2022 and 2023.

The noBB runs did NOT remove all fire emissions: they read the
`ne30np4_CONUSlandMasked_80kmBuffer` copies, in which only CONUS land + an
80 km buffer is zeroed. Fires outside that region (Canada, Alaska, Mexico)
are present in BOTH runs, so their O3 cancels in BASE - noBB. This script
shows what was in each run and what was actually removed.

Outputs (EMISSIONS_CHECK_FIGURE_DIR):
  BB_<SP>_meanflux_AprOct_2022-2023.png   per species, rows 2022 / 2023 /
        2022-2023, columns BASE input / noBB input / removed (BASE - noBB)
  BB_daily_totals_by_region_AprOct.png    daily emission totals, removed region
        vs outside it north of 49N (Canada + Alaska) vs the rest of the box
  BB_emission_totals_AprOct.csv           Apr-Oct totals (Gg) per species/year/region

Mean flux is plotted on the native ne30np4 cells (SCRIP polygons), in
kg km-2 day-1. Run in the `base` env (matplotlib + cartopy).

    python plot_bb_emissions_AprOct_BGO3.py                 # CO NO NO2 CH2O
    python plot_bb_emissions_AprOct_BGO3.py --species CO NO

    8 Oct 2026: VERSION 1.0
"""
import argparse, sys, pathlib
import numpy as np, pandas as pd, xarray as xr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm
import cartopy.crs as ccrs, cartopy.feature as cfeat

_ROOT = next(p for p in pathlib.Path(__file__).resolve().parents
             if (p / 'config' / 'paths.py').exists())
sys.path.insert(0, str(_ROOT))
import config  # noqa: F401
from config.paths import (BB_EMIS_NE30NP4_DIR, BB_EMIS_NE30NP4_CONUSMASKED_DIR,
                          SCRIP_NE30NP4, MASK_NE30NP4_CONUS_80KM,
                          EMISSIONS_CHECK_FIGURE_DIR, BGO3_YEARS,
                          BGO3_START_MMDD, BGO3_END_MMDD, ensure_dir)

FNAME = "qfed.emis_{sp}_bb_surface_daily_20171201T20231231_ne30np4_mol_c20240126.nc"
EXT = [-150, -55, 18, 68]            # North America: shows Canadian fires
R_EARTH = 6.37122e6                  # m, CESM shr_const_rearth
NA = 6.02214076e23
# molecules cm-2 s-1 -> kg km-2 day-1 is MW[g/mol]/NA * 1e10 cm2/km2 * 86400 s / 1e3
TO_KG_KM2_DAY = 1e10 * 86400 / 1e3 / NA
# g/mol. Not taken from the files: the NO and NO2 files carry a global
# molecular_weight of 184 and no emiss attribute.
MW = {"CO": 28.01, "NO": 30.01, "NO2": 46.01, "CH2O": 30.03, "SO2": 64.07,
      "NH3": 17.03, "C2H6": 30.07, "C3H8": 44.10, "CH3OH": 32.04, "ISOP": 68.12}

ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
ap.add_argument("--species", nargs="+", default=["CO", "NO", "NO2", "CH2O"])
args = ap.parse_args()
FIG = ensure_dir(EMISSIONS_CHECK_FIGURE_DIR)

# --- grid, cell area, zeroed-region mask ------------------------------------
g = xr.open_dataset(SCRIP_NE30NP4)
clon = g.grid_corner_lon.values.copy(); clat = g.grid_corner_lat.values
clon = np.where(clon > 180, clon - 360, clon)
lon = g.grid_center_lon.values; lon = np.where(lon > 180, lon - 360, lon)
lat = g.grid_center_lat.values
area_km2 = g.grid_area.values * R_EARTH**2 / 1e6          # radians^2 -> km2
inbox = ((lon >= EXT[0] - 5) & (lon <= EXT[1] + 5) & (lat >= EXT[2] - 5) &
         (lat <= EXT[3] + 5) & (np.ptp(clon, axis=1) < 180))
polys = np.stack([clon[inbox], clat[inbox]], axis=-1)

mk = xr.open_dataset(MASK_NE30NP4_CONUS_80KM)
zeroed = mk[list(mk.data_vars)[0]].values.astype(bool)     # True = zeroed in noBB
keep = ~zeroed
NA_BOX = (lon >= -170) & (lon <= -50) & (lat >= 10) & (lat <= 75)
REGIONS = {"removed (CONUS + 80 km)": zeroed,
           "kept, N of 49N (Canada+Alaska)": NA_BOX & keep & (lat >= 49),
           "kept, rest of N. America box": NA_BOX & keep & (lat < 49)}


def season_idx(date, yr):
    lo = int(f"{yr}{BGO3_START_MMDD.replace('-', '')}")
    hi = int(f"{yr}{BGO3_END_MMDD.replace('-', '')}")
    return np.where((date >= lo) & (date <= hi))[0]


def read(d, sp):
    """Daily flux (kg km-2 day-1) for each ozone season; dict yr -> (dates, arr)."""
    ds = xr.open_dataset(d / FNAME.format(sp=sp), decode_times=False)
    mw = MW.get(sp) or float(ds.emiss.attrs.get("molecular_weight",
                                               ds.attrs["molecular_weight"]))
    date = ds.date.values
    out = {}
    for yr in BGO3_YEARS:
        i = season_idx(date, yr)
        out[yr] = (date[i], ds.emiss.isel(time=i).values * mw * TO_KG_KM2_DAY)
    return out


def deco(ax, title):
    ax.add_feature(cfeat.COASTLINE, lw=.4); ax.add_feature(cfeat.BORDERS, lw=.4)
    ax.add_feature(cfeat.STATES, lw=.2, edgecolor="gray")
    ax.tricontour(lon[inbox], lat[inbox], zeroed[inbox].astype(float), levels=[.5],
                  colors="k", linewidths=.8, linestyles="--", transform=ccrs.PlateCarree())
    ax.set_extent(EXT, ccrs.PlateCarree()); ax.set_title(title, fontsize=10)


def cellmap(ax, field, norm, cmap):
    pc = PolyCollection(polys, array=np.ma.masked_less_equal(field[inbox], 0),
                        cmap=cmap, norm=norm, edgecolors="none",
                        transform=ccrs.PlateCarree())
    ax.add_collection(pc); return pc


rows, daily = [], {}
for sp in args.species:
    full = read(BB_EMIS_NE30NP4_DIR, sp)
    masked = read(BB_EMIS_NE30NP4_CONUSMASKED_DIR, sp)
    for yr in BGO3_YEARS:                               # mask must match the file
        f, m = full[yr][1], masked[yr][1]
        left = np.nansum(m[:, zeroed]) / max(np.nansum(f[:, zeroed]), 1e-30)
        assert left < 1e-6, f"{sp} {yr}: masked file keeps {left:.2%} inside mask"
        assert np.allclose(f[:, keep], m[:, keep], equal_nan=True), f"{sp} {yr}: differ outside mask"
        for reg, sel in REGIONS.items():
            tot = (f[:, sel] * area_km2[sel]).sum(1) / 1e6      # Gg day-1
            daily[(sp, yr, reg)] = pd.Series(tot, index=pd.to_datetime(full[yr][0].astype(str)))
            rows.append(dict(species=sp, year=yr, region=reg, total_Gg=tot.sum()))
        rows.append(dict(species=sp, year=yr, region="N. America box total",
                         total_Gg=(f[:, NA_BOX] * area_km2[NA_BOX]).sum() / 1e6))

    means = {str(yr): (full[yr][1].mean(0), masked[yr][1].mean(0)) for yr in BGO3_YEARS}
    cat = lambda k: np.concatenate([full[y][1] if k == 0 else masked[y][1]
                                    for y in BGO3_YEARS]).mean(0)
    means[f"{BGO3_YEARS[0]}-{BGO3_YEARS[-1]}"] = (cat(0), cat(1))
    vmax = np.nanpercentile(means[f"{BGO3_YEARS[0]}-{BGO3_YEARS[-1]}"][0][inbox], 99.9)
    norm = LogNorm(vmin=vmax / 1e3, vmax=vmax)
    cmap = plt.get_cmap("YlOrRd").copy(); cmap.set_bad("white")

    fig, axs = plt.subplots(len(means), 3, figsize=(15, 3.3 * len(means)),
                            subplot_kw={"projection": ccrs.PlateCarree()},
                            constrained_layout=True)
    for r, (lab, (b, n)) in enumerate(means.items()):
        for c, (fld, ttl) in enumerate([(b, "BASE input (all fires)"),
                                        (n, "noBB input (CONUS+80 km zeroed)"),
                                        (b - n, "Removed = BASE $-$ noBB")]):
            pc = cellmap(axs[r, c], fld, norm, cmap)
            deco(axs[r, c], f"{ttl}\nApr$-$Oct {lab}")
    fig.colorbar(pc, ax=axs, orientation="horizontal", shrink=.5, pad=.02,
                 extend="both").set_label(f"Mean BB {sp} emission (kg km$^{{-2}}$ day$^{{-1}}$)")
    fig.suptitle(f"QFED2.6 biomass-burning {sp} emissions used by CONUSBGO3 "
                 f"(dashed = zeroed region in noBB)", fontsize=13, fontweight="bold")
    out = FIG / f"BB_{sp}_meanflux_AprOct_{BGO3_YEARS[0]}-{BGO3_YEARS[-1]}.png"
    fig.savefig(out, dpi=150); plt.close(fig); print("saved", out)

# --- daily totals by region ---------------------------------------------------
fig, axs = plt.subplots(len(args.species), len(BGO3_YEARS), sharex="col",
                        figsize=(6.5 * len(BGO3_YEARS), 2.6 * len(args.species)),
                        squeeze=False, constrained_layout=True)
for i, sp in enumerate(args.species):
    for j, yr in enumerate(BGO3_YEARS):
        ax = axs[i, j]
        for reg, col in zip(REGIONS, ["C3", "C0", "0.5"]):
            s = daily[(sp, yr, reg)]
            ax.plot(s.index, s.values, color=col, lw=1, label=f"{reg}: {s.sum():.0f} Gg")
        ax.set_title(f"{sp} {yr}", fontsize=10); ax.set_ylabel("Gg day$^{-1}$")
        ax.legend(fontsize=7, frameon=False, loc="upper left")
fig.suptitle("Daily QFED2.6 biomass-burning emissions, Apr$-$Oct (red = removed in noBB)",
             fontsize=12, fontweight="bold")
out = FIG / "BB_daily_totals_by_region_AprOct.png"
fig.savefig(out, dpi=150); plt.close(fig); print("saved", out)

df = pd.DataFrame(rows)
out = FIG / "BB_emission_totals_AprOct.csv"
df.to_csv(out, index=False, float_format="%.2f"); print("saved", out)
print(df.pivot_table(index=["species", "region"], columns="year", values="total_Gg").round(1))
