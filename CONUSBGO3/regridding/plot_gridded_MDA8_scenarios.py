#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_gridded_MDA8_scenarios.py

Figure of the unified 1x1 gridded MDA8 O3 deliverable (produced by
Regrid_ne30_surfO3_to_1x1_conserve.py): a 3-panel CONUS map of the seasonal-mean
MDA8 O3 over Apr-Oct 2022-2023 —
  (a) BASE (absolute), (b) BASE - noAnthro, (c) BASE - noBB (differences).
Panels (b)/(c) use independent diverging scales (anthro effect ~= 8x the BB effect).

Run in the `base` env (matplotlib + cartopy). Edit MDA8_FILE / period as needed.

    Madankui Tao, 8 Jul 2026
"""
import numpy as np, xarray as xr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import cartopy.crs as ccrs, cartopy.feature as cfeat

OUT="/net/fs09/d0/taoma528/ProcessedData/DanJaffeMUSICAPostprocessing/"
MDA8_FILE=OUT+"Regridded1deg/MUSICAv0_ne30_CONUS1x1_MDA8O3_BGO3_2022-2023_AprOct_MTao_c20260708.nc"
FIG="/net/fs09/d0/taoma528/Figures/CESM_analysis/BGO3/regrid_trial/"
EXT=[-125,-66,23,50]

m=xr.open_dataset(MDA8_FILE); M=m["MDA8O3"]; lon=m.lon.values; lat=m.lat.values
si={str(s):i for i,s in enumerate(m.scenario.values)}
BASE=M.isel(scenario=si["BASE"]).mean("time")           # 2-yr Apr-Oct mean
dA=BASE-M.isel(scenario=si["noAnthro"]).mean("time")    # anthropogenic contribution
dB=BASE-M.isel(scenario=si["noBB"]).mean("time")        # biomass-burning contribution
vA=float(np.ceil(np.nanpercentile(np.abs(dA),99)))
vB=float(np.ceil(np.nanpercentile(np.abs(dB),99)))

def deco(ax,title):
    ax.add_feature(cfeat.COASTLINE,lw=.5); ax.add_feature(cfeat.STATES,lw=.3,edgecolor="gray")
    ax.add_feature(cfeat.BORDERS,lw=.4); ax.set_extent(EXT,ccrs.PlateCarree())
    ax.set_title(title,fontsize=12,fontweight="bold")

fig,axs=plt.subplots(1,3,figsize=(15,3.7),subplot_kw={"projection":ccrs.PlateCarree()},
                     constrained_layout=True)
p0=axs[0].pcolormesh(lon,lat,BASE,cmap="YlOrRd",vmin=25,vmax=60,shading="auto",transform=ccrs.PlateCarree())
deco(axs[0],"(a) BASE")
fig.colorbar(p0,ax=axs[0],orientation="horizontal",pad=0.03,aspect=28,extend="both").set_label("MDA8 O$_3$ (ppb)")
p1=axs[1].pcolormesh(lon,lat,dA,cmap="RdBu_r",vmin=-vA,vmax=vA,shading="auto",transform=ccrs.PlateCarree())
deco(axs[1],"(b) BASE $-$ noAnthro")
fig.colorbar(p1,ax=axs[1],orientation="horizontal",pad=0.03,aspect=28,extend="both").set_label("$\\Delta$MDA8 O$_3$ (ppb)")
p2=axs[2].pcolormesh(lon,lat,dB,cmap="RdBu_r",vmin=-vB,vmax=vB,shading="auto",transform=ccrs.PlateCarree())
deco(axs[2],"(c) BASE $-$ noBB")
fig.colorbar(p2,ax=axs[2],orientation="horizontal",pad=0.03,aspect=28,extend="both").set_label("$\\Delta$MDA8 O$_3$ (ppb)")
fig.suptitle("MUSICAv0 ne30 (1x1, mass-conservative) seasonal-mean MDA8 O$_3$, Apr-Oct 2022-2023",
             fontsize=13,fontweight="bold")
out=FIG+"MDA8_BASE_and_diffs_seasonmean_2022-2023.png"
fig.savefig(out,dpi=140); print("saved",out)
