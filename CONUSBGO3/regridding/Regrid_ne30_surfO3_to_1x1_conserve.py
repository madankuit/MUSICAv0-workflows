#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Regrid_ne30_surfO3_to_1x1_conserve.py

Mass-conservative (ESMF first-order) regrid of MUSICAv0 ne30np4 surface O3 to a
regular 1x1 deg CONUS grid, for the CONUS Background-O3 (Dan Jaffe) experiments
(BASE / noAnthro / noBB, 2022 & 2023, Apr-Oct). Produces:

  Regridded1deg/hourly/CONUS1x1_UTChourlySurfO3.<label>.<start>T<end>.nc   (6 files, ppb)
  Regridded1deg/MUSICAv0_ne30_CONUS1x1_MDA8O3_BGO3_2022-2023_AprOct_MTao_c<YMD>.nc
        -> MDA8O3(scenario, time, lat, lon)  [ppb]  (unified, shareable)

Conservative weights (ne30np4 -> 1x1) are read from CESM22/grids/; they were
generated once with esmpy 8.7 (rootxesmf env). Application here is a sparse
mat-mul (base env). MDA8 mirrors the point pipeline
(Extract_givenmonitorO3_hourly_dailyMDA8_toNetCDF.py): EPA 8-h rolling means,
per-grid-cell summertime UTC offset via timezonefinder.

MODIFICATION HISTORY:
    Madankui Tao, 8 Jul 2026: VERSION 1.0
"""
import os, datetime
import numpy as np, pandas as pd, xarray as xr
from scipy.sparse import coo_matrix
from timezonefinder import TimezoneFinder
from datetime import datetime as dtmod
import pytz

# ================= paths =================
G    = "/net/fs09/d0/taoma528/CESM22/grids/"
OUT  = "/net/fs09/d0/taoma528/ProcessedData/DanJaffeMUSICAPostprocessing/"
ARCH = "/net/fs09/d0/taoma528/CESM22/archive/"
DST  = OUT + "Regridded1deg/"
HRLY = DST + "hourly/"
os.makedirs(HRLY, exist_ok=True)
WGT  = G + "ESMFmap_ne30np4_TO_1x1_conserve_c20260708.nc"
h2ex = ARCH + ("f.e22.FCnudged.ne30_ne30_mg17.BGO3.BASEY20220401TY20230401/atm/hist/"
               "f.e22.FCnudged.ne30_ne30_mg17.BGO3.BASEY20220401TY20230401.cam.h2.2022-05-04-03600.nc")

# ================= cases =================
CASE = {
 ('BASE',2022):    'f.e22.FCnudged.ne30_ne30_mg17.BGO3.BASEY20220401TY20230401',
 ('BASE',2023):    'f.e22.FCnudged.ne30_ne30_mg17.BGO3.BASEY20230401TY20231101',
 ('noAnthro',2022):'f.e22.FCnudged.ne30_ne30_mg17.BGO3.noANTHROemisCONUS80kmBufferY20220401TY20221101',
 ('noAnthro',2023):'f.e22.FCnudged.ne30_ne30_mg17.BGO3.noANTHROemisCONUS80kmBufferY20230401TY20231101',
 ('noBB',2022):    'f.e22.FCnudged.ne30_ne30_mg17.BGO3.noBBemisCONUS80kmBufferY20220401TY20221101',
 ('noBB',2023):    'f.e22.FCnudged.ne30_ne30_mg17.BGO3.noBBemisCONUS80kmBufferY20230401TY20231101',
}
LABEL = {('BASE',2022):'BASE2022',('BASE',2023):'BASE2023',
         ('noAnthro',2022):'noAnthro2022',('noAnthro',2023):'noAnthro2023',
         ('noBB',2022):'noBB2022',('noBB',2023):'noBB2023'}
def merged_path(cn):
    import glob
    hits = glob.glob(OUT+f"h2_surfO3_merged/{cn}.cam.h2.surflev.O3.*.nc")
    if not hits: raise FileNotFoundError(cn)
    return sorted(hits)[-1]

SCEN = ['BASE','noAnthro','noBB']; YEARS = [2022,2023]
STARTMMDD, ENDMMDD = '04-01', '10-31'

# ================= CONUS 1x1 target =================
LATMIN,LATMAX,LONMIN,LONMAX = 24,50,235,294   # CONUS; lon 0..360 (=-125..-66)
latc = np.arange(LATMIN,LATMAX+1,1); lonc = np.arange(LONMIN,LONMAX+1,1)
nlat,nlon = latc.size, lonc.size

# ================= build conservative operator restricted to CONUS =================
# Global dest grid (FV1x1grid_info): lat -90..90 (181), lon 0..360 (361).
lat_grid=np.arange(-90,91,1.0); lon_grid=np.arange(0,361,1.0)
NLATG,NLONG=lat_grid.size,lon_grid.size; NGLOB=NLATG*NLONG; NSRC=48602
w=xr.open_dataset(WGT); col=w['col'].values-1; row=w['row'].values-1; S=w['S'].values
Mg=coo_matrix((S,(row,col)),shape=(NGLOB,NSRC)).tocsr()
wsum=np.asarray(Mg.sum(1)).ravel()
slat=xr.open_dataset(h2ex)['lat'].values; slon=xr.open_dataset(h2ex)['lon'].values
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
sel=np.where((wsum>0.999)&(tlat>=LATMIN)&(tlat<=LATMAX)&(tlon>=LONMIN)&(tlon<=LONMAX))[0]
iy=(tlat[sel]-LATMIN).astype(int); ix=(tlon[sel]-LONMIN).astype(int)   # exact bijection
Mc=Mg[sel,:]; wsum_c=wsum[sel]; ncell=sel.size
cell_lat=latc[iy]; cell_lon180=np.where(lonc[ix]>180,lonc[ix]-360,lonc[ix])
assert len(set(zip(iy.tolist(),ix.tolist())))==ncell, "placement not unique!"
print(f"CONUS cells mapped: {ncell} (grid {nlat}x{nlon}={nlat*nlon})")

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
YMD = dtmod.now().strftime('%Y%m%d')
def prov(extra):
    a = dict(
        title="MUSICAv0 ne30 CONUS Background-O3: surface O3 regridded to 1x1 deg (mass-conservative)",
        project="CONUS Background Ozone (Dan Jaffe collaboration)",
        source=("MUSICAv0 = CESM2.2 CAM-chem (MOZART TS1), ne30np4 (~111 km global) "
                "spectral-element, FCnudged nudged to MERRA-2"),
        institution="Smithsonian Astrophysical Observatory (Center for Astrophysics | Harvard & Smithsonian)",
        scenarios=("BASE = all emissions; noAnthro = CONUS anthropogenic emissions zeroed (land, 80 km buffer); "
                   "noBB = CONUS biomass-burning emissions zeroed (land, 80 km buffer)"),
        case_names="; ".join(f"{LABEL[k]}: {v}" for k,v in CASE.items()),
        horizontal_regrid="ESMF first-order conservative remap, ne30np4 -> 1x1 deg (dest weight-sums = 1)",
        regrid_weight_file=WGT,
        regrid_source_scrip=G+"ne30np4_091226_pentagons.nc",
        regrid_dest_grid=G+"FV1x1grid_info_c20241105.nc",
        processed_by="Madankui Tao (taoma528)",
        processing_date=dtmod.now().strftime('%Y-%m-%d %H:%M:%S'),
        machine="MIT Svante (svante9.mit.edu)",
        conda_env="base (application/MDA8); weights generated with esmpy 8.7 in env rootxesmf",
        contact="taoma528@mit.edu",
        Conventions="CF-1.9",
    )
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
lon_out=np.where(lonc>180,lonc-360.0,lonc).astype('float64')     # -128..-64 for output
mda8_all={}   # (scenario,year) -> (214, ncell)
mda8_dates={}
for sc in SCEN:
    for yr in YEARS:
        cn=CASE[(sc,yr)]; lab=LABEL[(sc,yr)]; p=merged_path(cn)
        print(f"\n=== {lab} === {os.path.basename(p)}")
        times,hourly=regrid_hourly(p)                            # (nt, ncell) ppb, UTC

        # ---- save hourly gridded (season subset 04-01..10-31 UTC) ----
        tmask=(pd.DatetimeIndex(times)>=pd.Timestamp(f'{yr}-{STARTMMDD}T00')) & \
              (pd.DatetimeIndex(times)<=pd.Timestamp(f'{yr}-{ENDMMDD}T23'))
        hmap=to_map(hourly[tmask]); htimes=times[tmask]
        dsh=xr.Dataset(
            {'O3':(('time','lat','lon'),hmap,
                   {'units':'ppb','long_name':'Hourly surface O3 (UTC), conservatively regridded to 1x1'})},
            coords={'time':htimes,'lat':latc.astype('float64'),'lon':lon_out},
            attrs=prov({'title':f'MUSICAv0 ne30 {lab} hourly surface O3 on 1x1 CONUS grid (UTC)',
                        'case_name':cn,'scenario':sc,'year':yr,
                        'time_note':'UTC timestamps; hours over Apr 1 - Oct 31'}))
        dsh['lat'].attrs.update(units='degrees_north',standard_name='latitude')
        dsh['lon'].attrs.update(units='degrees_east',standard_name='longitude')
        fpath=HRLY+f'CONUS1x1_UTChourlySurfO3.{lab}.{yr}-{STARTMMDD}T{yr}-{ENDMMDD}.MTao_c{YMD}.nc'
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

# ================= unified MDA8 file =================
time_all=pd.DatetimeIndex(np.concatenate([mda8_dates[y].values for y in YEARS]))
arr=np.full((len(SCEN),len(time_all),nlat,nlon),np.nan,dtype='float32')
for si,sc in enumerate(SCEN):
    blk=0
    for yr in YEARS:
        n=len(mda8_dates[yr])
        arr[si,blk:blk+n]=to_map(mda8_all[(sc,yr)]); blk+=n
dsm=xr.Dataset(
    {'MDA8O3':(('scenario','time','lat','lon'),arr,
        {'units':'ppb','long_name':'Daily maximum 8-hour average surface O3 (MDA8), local-time dates',
         'cell_methods':'time: mean (interval: 8 hours) time: maximum within days'})},
    coords={'scenario':np.array(SCEN),'time':time_all,
            'lat':latc.astype('float64'),'lon':lon_out},
    attrs=prov({'title':'MUSICAv0 ne30 CONUS Background-O3: 1x1 deg daily MDA8 O3, 3 scenarios, 2022-2023 (Apr-Oct)',
                'mda8_method':('EPA convention: 8-h rolling mean (>=6 valid hours); daily max over windows '
                               'ending 07-23 local; day valid if >=13 of 17 windows'),
                'local_time':'per-grid-cell summertime (DST) UTC offset from timezonefinder',
                'scenario_dim':'BASE, noAnthro, noBB',
                'time_note':'428 daily dates = Apr1-Oct31 2022 (214) + Apr1-Oct31 2023 (214)'}))
dsm['lat'].attrs.update(units='degrees_north',standard_name='latitude')
dsm['lon'].attrs.update(units='degrees_east',standard_name='longitude')
dsm['scenario'].attrs.update(long_name='emission scenario')
fpath=DST+f'MUSICAv0_ne30_CONUS1x1_MDA8O3_BGO3_2022-2023_AprOct_MTao_c{YMD}.nc'
dsm.to_netcdf(fpath,format='NETCDF4',
              encoding={'MDA8O3':{'zlib':True,'complevel':4,'_FillValue':np.float32(np.nan)}})
print("\nUNIFIED MDA8 ->",fpath)
print("done.")
