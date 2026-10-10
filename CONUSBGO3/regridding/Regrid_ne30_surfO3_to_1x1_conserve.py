#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Regrid_ne30_surfO3_to_1x1_conserve.py

Mass-conservative (ESMF first-order) regrid of MUSICAv0 ne30np4 surface O3 to a
regular 1x1 deg CONUS grid, for the CONUS background-O3 experiments
(BASE / noAnthro / noBBCONUS / noBBGlobal, 2022 & 2023, Apr-Oct). Produces:

  Regridded1deg/hourly/CONUS1x1_UTChourlySurfO3.<label>.<start>T<end>.<tag>_c<YMD>.nc  (ppb)
  Regridded1deg/MUSICAv0_ne30_CONUS1x1_MDA8O3_2022-2023_AprOct_c<YMD>.nc
        -> MDA8O3(scenario, time, lat, lon)  [ppb]  (unified, shareable)

With --append-to, only --scenarios are computed and the other scenarios are copied
unchanged from an existing unified file (its 'noBB' is renamed noBBCONUS).

Conservative weights (ne30np4 -> 1x1) are read from CESM22/grids/; they were
generated once with esmpy 8.7 (rootxesmf env). Application here is a sparse
mat-mul (base env). MDA8 mirrors the point pipeline
(Extract_givenmonitorO3_hourly_dailyMDA8_toNetCDF.py): EPA 8-h rolling means,
per-grid-cell summertime UTC offset via timezonefinder.

All paths, case names and the provenance tag come from config/paths.py.

MODIFICATION HISTORY:
    8 Jul 2026: VERSION 1.0
    31 Aug 2026: VERSION 1.1
    - Paths, case names and provenance moved to config/paths.py
    9 Oct 2026: VERSION 1.2
    - 4th scenario noBBGlobal; 'noBB' renamed noBBCONUS. Merged files are found
      per (case, year) because the noBBGlobal case serves both years
    - --scenarios / --append-to / --cdate; --check-against for a regression test

    # add noBBGlobal to the delivered 3-scenario file
    python Regrid_ne30_surfO3_to_1x1_conserve.py --scenarios noBBGlobal \
        --append-to <old unified file> --cdate 20261010
    # regression: recompute BASE 2022 and require equality with the old file
    python Regrid_ne30_surfO3_to_1x1_conserve.py --scenarios BASE --years 2022 \
        --check-against <old unified file> --no-write
    - --experiment ne0CONUS (ne0CONUSne30x8 runs) and --target 1x1 | 0.15: the SE
      grid, weights, destination grid and output names come from
      BGO3_EXPERIMENTS; the MDA8 code is unchanged

    python Regrid_ne30_surfO3_to_1x1_conserve.py --experiment ne0CONUS --target 1x1
    python Regrid_ne30_surfO3_to_1x1_conserve.py --experiment ne0CONUS --target 0.15
"""
import os, datetime
import numpy as np, pandas as pd, xarray as xr
from scipy.sparse import coo_matrix
from timezonefinder import TimezoneFinder
from datetime import datetime as dtmod
import pytz

# ================= configuration =================
# Every path, case name and provenance string comes from config/paths.py.
import sys, pathlib, glob
_ROOT = next(p for p in pathlib.Path(__file__).resolve().parents
             if (p / 'config' / 'paths.py').exists())
sys.path.insert(0, str(_ROOT))
import config  # noqa: F401  - also puts functions/ on sys.path
from config.paths import (
    BGO3_EXPERIMENTS, BGO3_START_MMDD, BGO3_END_MMDD,
    AUTHOR_TAG, PROCESSED_BY, CONTACT, INSTITUTION, MACHINE,
    bgo3_merged_surfo3_path, case_hist_dir, ensure_dir,
)
import argparse
ap = argparse.ArgumentParser()
ap.add_argument('--experiment', default='ne30', choices=list(BGO3_EXPERIMENTS))
ap.add_argument('--target', default='1x1', help='destination grid key in BGO3_EXPERIMENTS[...]["targets"]')
ap.add_argument('--scenarios', nargs='+', default=None,
                help='scenarios to compute (hourly files are written for these)')
ap.add_argument('--years', nargs='+', type=int, default=None)
ap.add_argument('--append-to', default=None,
                help='existing unified MDA8 file; scenarios not computed are copied from it')
ap.add_argument('--check-against', default=None,
                help='unified MDA8 file; computed MDA8 must equal it (regression test)')
ap.add_argument('--no-write', action='store_true', help='write nothing (use with --check-against)')
ap.add_argument('--cdate', default=None, help='YYYYMMDD stamp for output names (default today)')
args = ap.parse_args()
OLD_NAMES = {'noBB': 'noBBCONUS'}         # scenario names in files delivered before Oct 2026
EXP = BGO3_EXPERIMENTS[args.experiment]
assert args.target in EXP['targets'], f"--target must be one of {list(EXP['targets'])}"
GRIDINFO, WGT, GTAG = EXP['targets'][args.target]          # e.g. FV1x1 gridinfo, weights, 'CONUS1x1'
GRIDINFO, WGT = str(GRIDINFO), str(WGT)

DST  = str(ensure_dir(EXP['regrid_dir'])) + "/"
HRLY = str(ensure_dir(EXP['hourly_dir'])) + "/"

# ================= cases =================
CASE  = dict(EXP['cases'])
LABEL = {k: f"{EXP['label_prefix']}{k[0]}{k[1]}" for k in CASE}   # e.g. noBBGlobal2023

# Any h2 file of the first BASE case provides the SE coordinates / grid metadata.
_h2ex_hist = case_hist_dir(CASE[('BASE', EXP['years'][0])])
_h2ex_hits = sorted(glob.glob(str(_h2ex_hist / '*.cam.h2.*.nc')))
if not _h2ex_hits:
    raise FileNotFoundError(f'No h2 history files found under {_h2ex_hist}')
h2ex = _h2ex_hits[0]

SCEN = list(EXP['scenarios']); YEARS = list(EXP['years'])
args.scenarios = args.scenarios or SCEN; args.years = args.years or YEARS
assert set(args.scenarios) <= set(SCEN) and set(args.years) <= set(YEARS)
RUN_SCEN = [sc for sc in SCEN if sc in args.scenarios]; RUN_YEARS = [y for y in YEARS if y in args.years]
STARTMMDD, ENDMMDD = BGO3_START_MMDD, BGO3_END_MMDD

# ================= CONUS box on the destination grid =================
LATMIN,LATMAX,LONMIN,LONMAX = 24,50,-125,-66   # CONUS; cell centres inside, lon -180..180

# ================= build conservative operator restricted to CONUS =================
# Global dest grid from its gridinfo file (1x1: lat -90..90 (181), lon 0..360 (361);
# 0.15: lat -89.95..89.9 (1200), lon -179.95..179.9 (2400)).
_gi=xr.open_dataset(GRIDINFO)
lat_grid=_gi['lat'].values.astype('float64'); lon_grid=_gi['lon'].values.astype('float64')
NLATG,NLONG=lat_grid.size,lon_grid.size; NGLOB=NLATG*NLONG; NSRC=EXP['ncol']
w=xr.open_dataset(WGT); col=w['col'].values-1; row=w['row'].values-1; S=w['S'].values
assert col.max()+1==NSRC and row.max()+1<=NGLOB, f"{WGT}: does not map {EXP['grid']} -> {GRIDINFO}"
Mg=coo_matrix((S,(row,col)),shape=(NGLOB,NSRC)).tocsr()
wsum=np.asarray(Mg.sum(1)).ravel()
slat=xr.open_dataset(h2ex)['lat'].values; slon=xr.open_dataset(h2ex)['lon'].values
assert slat.size==NSRC, (h2ex, slat.size)
clat=Mg.dot(slat); clon=Mg.dot(slon)                      # area-wtd centroid (for ORDER detection only)

# --- detect ESMF flattened dest ordering (C: lon-fastest vs F: lat-fastest) ---
good=np.where(wsum>0.5)[0]
def _res(il,io):
    dlon=np.abs(((lon_grid[io]-clon[good]+180)%360)-180)
    return np.abs(lat_grid[il]-clat[good]).mean()+dlon.mean()
rC=_res(good//NLONG, good%NLONG); rF=_res(good%NLATG, good//NLATG)
corder=rC<=rF
print(f"dest ordering: C(lon-fast) res={rC:.3f}  F(lat-fast) res={rF:.3f} -> use {'C' if corder else 'F'}")
allr=np.arange(NGLOB)
ILAT=(allr//NLONG) if corder else (allr%NLATG)
ILON=(allr%NLONG)  if corder else (allr//NLATG)
tlat=lat_grid[ILAT]; tlon=lon_grid[ILON]                  # TRUE center of each global dest cell
tlon180=np.where(tlon>180,tlon-360,tlon)
sel=np.where((wsum>0.999)&(tlat>=LATMIN)&(tlat<=LATMAX)&(tlon180>=LONMIN)&(tlon180<=LONMAX))[0]
latc=np.unique(tlat[sel]); lonc=np.unique(tlon180[sel]); nlat,nlon=latc.size,lonc.size
iy=np.searchsorted(latc,tlat[sel]); ix=np.searchsorted(lonc,tlon180[sel])   # exact bijection
Mc=Mg[sel,:]; wsum_c=wsum[sel]; ncell=sel.size
cell_lat=latc[iy]; cell_lon180=lonc[ix]
assert len(set(zip(iy.tolist(),ix.tolist())))==ncell, "placement not unique!"
print(f"CONUS cells mapped: {ncell} (grid {nlat}x{nlon}={nlat*nlon}) "
      f"lat {latc[0]:.3f}..{latc[-1]:.3f} lon {lonc[0]:.3f}..{lonc[-1]:.3f}")

# ================= per-cell summertime UTC offset (timezonefinder) =================
tf=TimezoneFinder()
def summer_offset(lat,lon180):
    tz=tf.timezone_at(lat=lat,lng=lon180)
    if tz is None: return int(round(lon180/15.0))
    dt=pytz.timezone(tz).localize(dtmod(2022,7,1),is_dst=True)
    return int(dt.utcoffset().total_seconds()//3600)
offset=np.array([summer_offset(la,lo) for la,lo in zip(cell_lat,cell_lon180)])
print("cell UTC offsets present:",sorted(set(offset.tolist())))

# ================= MDA8 (EPA), vectorized over cells, per point pipeline =================
def compute_mda8(o3_hourly_utc, datesv1, utc_off, min_hours=6, min_blocks=13):
    t_local = pd.DatetimeIndex(o3_hourly_utc.time.values) + pd.to_timedelta(utc_off,'h')
    O3 = o3_hourly_utc.assign_coords(time=('time',t_local))
    O3_8h = O3.rolling(time=8, min_periods=min_hours).mean()
    O3_8h = O3_8h.where(O3_8h['time'].dt.hour.isin(np.arange(7,24)))
    day = pd.DatetimeIndex(O3_8h.time.values).floor('D')
    O3_8h = O3_8h.assign_coords(day=('time',day))
    mda8 = O3_8h.groupby('day').max('time',skipna=True)
    cnt  = O3_8h.groupby('day').count('time')
    mda8 = mda8.where(cnt>=min_blocks)
    tgt = pd.to_datetime(datesv1)
    return mda8.reindex(day=tgt).rename({'day':'time'}).assign_coords(time=('time',tgt))

# ================= common provenance attrs =================
YMD = args.cdate or dtmod.now().strftime('%Y%m%d')
GRES = {'1x1': '1x1 deg', '0.15': '0.15x0.15 deg'}.get(args.target, args.target)
SCEN_TEXT = {
    'ne30': None,   # filled below (kept verbatim from v1.2)
    'ne0CONUS': dict(
        scenarios=("BASE = all emissions (hourly NEI2022v2 over CONUS + CAMS-GLOB-ANT v6.2 elsewhere, "
                   "QFED2.6 hi-res biomass burning); noBBCONUS = QFED2.6 biomass-burning emissions zeroed over "
                   "CONUS land + 80 km buffer only (fires in Canada/Alaska/Mexico kept, so BASE - noBBCONUS = "
                   "O3 from US fires). No global no-BB run exists on this grid."),
        scenario_notes=("BASE and noBBCONUS (.01) are hybrid starts from the same state on 2023-01-01, so CONUS "
                        "fires have been off in noBBCONUS since Jan 1 (3-month spin-up before the Apr-Oct window). "
                        "noBBCONUS .01 branched to .02 at 2023-04-01 with an identical namelist; only .02 output "
                        "is used. BASE and noBBCONUS share grid, GEOS-5 nudging and all "
                        "non-fire emissions. Hourly values are means stamped at the END of each hour; every series "
                        "starts at Apr 1 01:00 UTC. This BASE is a different model setup from the ne30 CONUSBGO3 "
                        "BASE (inventories, resolution, nudging), so the two sets are not interchangeable."),
    ),
}
def prov(extra):
    st = SCEN_TEXT.get(args.experiment)
    a = dict(
        title=f"MUSICAv0 {EXP['grid']} CONUS Background-O3: surface O3 regridded to {GRES} (mass-conservative)",
        project="CONUS background ozone: emission-zeroing sensitivity experiments",
        source=EXP['model'],
        institution=INSTITUTION,
        scenarios=("BASE = all emissions; "
                   "noAnthro = CONUS anthropogenic emissions zeroed (land, 80 km buffer); "
                   "noBBCONUS = QFED2.6 biomass-burning emissions zeroed over CONUS land + 80 km buffer only "
                   "(fires in Canada/Alaska/Mexico kept, so BASE - noBBCONUS = O3 from US fires; "
                   "called noBB in files delivered before Oct 2026); "
                   "noBBGlobal = biomass-burning emissions zeroed globally (QFED2.6 plus CMIP6 fire DMS and "
                   "num_so4_a1), so BASE - noBBGlobal = O3 from all fires incl. transported smoke and the "
                   "hemispheric fire background"),
        scenario_notes=("noBBGlobal is one continuous run from 2022-04-01, so its 2023 season starts from a "
                        "no-fire state; BASE2023, noAnthro2023 and noBBCONUS2023 branch from BASE2022 at "
                        "2023-04-01. Grid, MERRA-2 nudging, all other emissions and the 2022 initial "
                        "condition are identical across scenarios."),
        case_names="; ".join(f"{LABEL[k]}: {v}" for k,v in CASE.items()),
        horizontal_regrid=f"ESMF first-order conservative remap, {EXP['grid']} -> {GRES} (dest weight-sums = 1)",
        regrid_weight_file=WGT,
        regrid_source_scrip=str(EXP['scrip']),
        regrid_dest_grid=GRIDINFO,
        processed_by=PROCESSED_BY,
        processing_date=dtmod.now().strftime('%Y-%m-%d %H:%M:%S'),
        machine=MACHINE,
        conda_env="base (application/MDA8); weights generated with esmpy 8.7 in env rootxesmf",
        contact=CONTACT,
        Conventions="CF-1.9",
    )
    if st: a.update(st)
    a.update(extra); return a

# ================= regrid one merged file (chunked mat-mul) =================
def regrid_hourly(path):
    ds=xr.open_dataset(path); o3=ds['O3']; nt=o3.sizes['time']; times=ds.time.values
    outc=np.empty((nt,ncell),dtype='float32'); step=744
    for s in range(0,nt,step):
        e=min(s+step,nt)
        chunk=o3.isel(time=slice(s,e)).values.astype('float32')   # (nc, ncol)
        d=Mc.dot(chunk.T)/wsum_c[:,None]                          # (ncell, nc)
        outc[s:e,:]=(d.T*1e9).astype('float32')                  # ppb
    return times, outc

def to_map(vec_cell_time):   # (ntime, ncell) -> (ntime, nlat, nlon)
    g=np.full((vec_cell_time.shape[0],nlat,nlon),np.nan,dtype='float32')
    g[:,iy,ix]=vec_cell_time
    return g

# ================= main =================
lon_out=np.round(lonc,4).astype('float64')                       # -125..-66 for output
latc=np.round(latc,4)
merged_used={}  # label -> merged input file
mda8_all={}   # (scenario,year) -> (214, ncell)
mda8_dates={}
for sc in RUN_SCEN:
    for yr in RUN_YEARS:
        cn=CASE[(sc,yr)]; lab=LABEL[(sc,yr)]; p=bgo3_merged_surfo3_path(cn,yr,EXP['merged_dir'])
        print(f"\n=== {lab} === {cn}\n    merged: {os.path.basename(p)}")
        merged_used[lab]=os.path.basename(p)
        times,hourly=regrid_hourly(p)                            # (nt, ncell) ppb, UTC
        _t=pd.DatetimeIndex(times)
        assert (_t.year==yr).all() and _t.min()<=pd.Timestamp(f'{yr}-{STARTMMDD}T01') \
            and _t.max()>=pd.Timestamp(f'{yr}-{ENDMMDD}T23'), f"{lab}: merged spans {_t.min()}..{_t.max()}"

        # ---- save hourly gridded (season subset 04-01..10-31 UTC) ----
        tmask=(pd.DatetimeIndex(times)>=pd.Timestamp(f'{yr}-{STARTMMDD}T00')) & \
              (pd.DatetimeIndex(times)<=pd.Timestamp(f'{yr}-{ENDMMDD}T23'))
        hmap=to_map(hourly[tmask]); htimes=times[tmask]
        dsh=xr.Dataset(
            {'O3':(('time','lat','lon'),hmap,
                   {'units':'ppb','long_name':f'Hourly surface O3 (UTC), conservatively regridded to {GRES}'})},
            coords={'time':htimes,'lat':latc.astype('float64'),'lon':lon_out},
            attrs=prov({'title':f"MUSICAv0 {EXP['grid']} {lab} hourly surface O3 on {GRES} CONUS grid (UTC)",
                        'case_name':cn,'scenario':sc,'year':yr,'merged_input_file':os.path.basename(p),
                        'time_note':'UTC timestamps; hours over Apr 1 - Oct 31'}))
        dsh['lat'].attrs.update(units='degrees_north',standard_name='latitude')
        dsh['lon'].attrs.update(units='degrees_east',standard_name='longitude')
        fpath=HRLY+f'{GTAG}_UTChourlySurfO3.{lab}.{yr}-{STARTMMDD}T{yr}-{ENDMMDD}.{AUTHOR_TAG}_c{YMD}.nc'
        if not args.no_write:
            assert not os.path.exists(fpath), f"{fpath} exists"
            dsh.to_netcdf(fpath,format='NETCDF4',
                          encoding={'O3':{'zlib':True,'complevel':4,'_FillValue':np.float32(np.nan)}})
            print("  hourly ->",fpath)

        # ---- MDA8 (use full hourly for local-time completeness; target 04-01..10-31) ----
        datesv1=pd.date_range(f'{yr}-{STARTMMDD}',f'{yr}-{ENDMMDD}',freq='D')
        da=xr.DataArray(hourly,dims=['time','cell'],coords={'time':times})
        mda8=np.full((len(datesv1),ncell),np.nan,dtype='float32')
        for off in sorted(set(offset.tolist())):
            idx=np.where(offset==off)[0]
            res=compute_mda8(da.isel(cell=idx),datesv1.strftime('%Y-%m-%d').tolist(),off)
            mda8[:,idx]=res.values
        mda8_all[(sc,yr)]=mda8; mda8_dates[yr]=datesv1
        print(f"  MDA8 ppb: min={np.nanmin(mda8):.1f} max={np.nanmax(mda8):.1f} nan%={100*np.isnan(mda8).mean():.0f}")

# ================= regression check against an existing unified file =================
def _open_unified(path):
    m=xr.open_dataset(path)
    return m.assign_coords(scenario=[OLD_NAMES.get(str(v),str(v)) for v in m.scenario.values])

if args.check_against:
    ref=_open_unified(args.check_against)
    for (sc,yr),mda8 in mda8_all.items():
        r=ref['MDA8O3'].sel(scenario=sc,time=str(yr)).values
        new=to_map(mda8)
        same=np.array_equal(np.isnan(r),np.isnan(new)) and np.array_equal(r[~np.isnan(r)],new[~np.isnan(new)])
        dmax=np.nanmax(np.abs(r-new))
        print(f"CHECK {LABEL[(sc,yr)]} vs {os.path.basename(args.check_against)}: "
              f"{'IDENTICAL' if same else 'DIFFERENT'} (max |diff| {dmax:.3g} ppb)")
        assert same, f"{LABEL[(sc,yr)]} differs from {args.check_against}"
if args.no_write:
    print("--no-write: done."); sys.exit(0)

# ================= unified MDA8 file =================
assert RUN_YEARS==YEARS, "the unified file needs both years"
time_all=pd.DatetimeIndex(np.concatenate([mda8_dates[y].values for y in YEARS]))
arr=np.full((len(SCEN),len(time_all),nlat,nlon),np.nan,dtype='float32')
copied=[]
if args.append_to:
    base=_open_unified(args.append_to)
    assert np.array_equal(pd.DatetimeIndex(base.time.values),time_all), "time axis differs from --append-to"
    assert np.array_equal(base.lat.values,latc) and np.array_equal(base.lon.values,lon_out), "grid differs"
for si,sc in enumerate(SCEN):
    if sc in RUN_SCEN:
        blk=0
        for yr in YEARS:
            n=len(mda8_dates[yr])
            arr[si,blk:blk+n]=to_map(mda8_all[(sc,yr)]); blk+=n
    else:
        assert args.append_to and sc in base.scenario.values, f"{sc}: not computed and not in --append-to"
        arr[si]=base['MDA8O3'].sel(scenario=sc).values; copied.append(sc)
case_names="; ".join(f"{LABEL[k]}: {v}" for k,v in CASE.items())
extra_src=(f"{', '.join(copied)} copied unchanged from {os.path.basename(args.append_to)} "
           f"(its 'noBB' renamed noBBCONUS); {', '.join(RUN_SCEN)} computed here" if copied else
           "all scenarios computed here")
dsm=xr.Dataset(
    {'MDA8O3':(('scenario','time','lat','lon'),arr,
        {'units':'ppb','long_name':'Daily maximum 8-hour average surface O3 (MDA8), local-time dates',
         'cell_methods':'time: mean (interval: 8 hours) time: maximum within days'})},
    coords={'scenario':np.array(SCEN),'time':time_all,
            'lat':latc.astype('float64'),'lon':lon_out},
    attrs=prov({'title':(f"MUSICAv0 {EXP['grid']} CONUS Background-O3: {GRES} daily MDA8 O3, {len(SCEN)} scenarios, "
                         f"{YEARS[0]}-{YEARS[-1]} (Apr-Oct)"),
                'case_names':case_names,
                'mda8_method':('EPA convention: 8-h rolling mean (>=6 valid hours); daily max over windows '
                               'ending 07-23 local; day valid if >=13 of 17 windows'),
                'local_time':'per-grid-cell summertime (DST) UTC offset from timezonefinder',
                'scenario_dim':', '.join(SCEN),
                'scenario_source':extra_src,
                'merged_input_files':'; '.join(f'{k}: {v}' for k,v in merged_used.items()),
                'time_note':(f"{len(time_all)} daily dates = " +
                             " + ".join(f"Apr1-Oct31 {y} ({len(mda8_dates[y])})" for y in YEARS))}))
dsm['lat'].attrs.update(units='degrees_north',standard_name='latitude')
dsm['lon'].attrs.update(units='degrees_east',standard_name='longitude')
dsm['scenario'].attrs.update(long_name='emission scenario')
fpath=DST+EXP['unified_stem'].format(tag=GTAG)+f'_c{YMD}.nc'
assert not os.path.exists(fpath), f"{fpath} exists; pass another --cdate or remove it"
dsm.to_netcdf(fpath,format='NETCDF4',
              encoding={'MDA8O3':{'zlib':True,'complevel':4,'_FillValue':np.float32(np.nan)}})
print("\nUNIFIED MDA8 ->",fpath)
print("done.")
