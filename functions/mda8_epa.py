#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mda8_epa.py

Daily maximum 8-hour average O3 (MDA8) following the 2015 O3 NAAQS data-handling
rules, 40 CFR Part 50 Appendix U, Section 3:

  (a) each hour is identified in LOCAL STANDARD TIME (LST, no daylight saving);
      hourly concentrations are truncated to 3 decimals in ppm (i.e. whole ppb)
  (b) 8-hour averages are stored at the START hour of the period; valid if >= 6 of
      the 8 hours are available (divisor 6 or 7). With <= 5 hours, zero is
      substituted for the missing hours and the average (divisor 8) is valid only
      if it exceeds the NAAQS level. Averages are truncated to 3 decimals in ppm.
  (c) the daily maximum is the highest of the 17 consecutive 8-hour averages
      beginning 07:00-15:00 through 23:00-07:00 (next day), LST
  (d) a daily maximum is valid if >= 13 of the 17 averages are valid, or if it
      is greater than the NAAQS level.

NAAQS level 0.070 ppm = 70 ppb; "greater than" a truncated whole-ppb value means
>= 71 ppb. Model output: CAM hourly means are stamped at the END of their hour, so
the hour-beginning label is stamp - 1 h. A tiny epsilon (1e-6 ppb) guards the
truncation against float round-off (e.g. 45.9999999 -> 46, not 45).

    from mda8_epa import mda8_epa, lst_offset_hours
    days, mda8 = mda8_epa(hourly_ppb, time_end_utc, lst_offsets, '2022-04-01', '2022-10-31')

    18 Oct 2026: VERSION 1.0
"""
from datetime import datetime
import numpy as np
import pandas as pd

NAAQS_PPB = 70.0
EPS = 1e-6


def lst_offset_hours(tz_name, lon, year=2023):
    """Standard-time UTC offset (h) of an IANA zone in `year` (no DST); nautical if no zone."""
    import pytz
    if tz_name is None:
        return int(round(lon / 15.0))
    loc = pytz.timezone(tz_name).localize(datetime(year, 7, 1), is_dst=False)
    return int((loc.utcoffset() - loc.dst()).total_seconds() // 3600)


def truncate_ppb(x):
    """Truncate to whole ppb (= 3 decimals in ppm), as Appendix U 3(a)/(b)."""
    return np.floor(np.asarray(x, dtype='float64') + EPS)


def mda8_epa(hourly_ppb, time_end_utc, lst_offsets, first_day, last_day, return_counts=False):
    """MDA8 (ppb, whole numbers; NaN = not valid) for local dates first_day..last_day.

    hourly_ppb   (ntime, ncell) hourly mean O3 in ppb
    time_end_utc (ntime,) datetime64, END of each hourly mean (CAM convention)
    lst_offsets  (ncell,) int hours, local standard time - UTC
    """
    H = np.asarray(hourly_ppb, dtype='float64')
    if H.ndim == 1:
        H = H[:, None]
    off = np.asarray(lst_offsets, dtype=int).reshape(-1)
    assert H.shape[1] == off.size, (H.shape, off.size)
    begin_utc = pd.DatetimeIndex(time_end_utc) - pd.Timedelta(hours=1)        # hour-beginning label
    days = pd.date_range(first_day, last_day)
    out = np.full((len(days), H.shape[1]), np.nan, dtype='float32')
    nval = np.zeros((len(days), H.shape[1]), dtype='int16')
    Ht = np.where(np.isfinite(H), truncate_ppb(H), np.nan)                     # 3(a) truncation
    for o in np.unique(off):
        idx = np.where(off == o)[0]
        lst = begin_utc + pd.Timedelta(hours=int(o))
        grid = pd.date_range(days[0] + pd.Timedelta(hours=7), days[-1] + pd.Timedelta(hours=30), freq='h')
        pos = grid.get_indexer(lst); ok = pos >= 0
        A = np.full((len(grid), idx.size), np.nan); A[pos[ok]] = Ht[ok][:, idx]
        v = np.nan_to_num(A); n = np.isfinite(A).astype(float)
        cs = np.vstack([np.zeros((1, idx.size)), np.cumsum(v, 0)]); cn = np.vstack([np.zeros((1, idx.size)), np.cumsum(n, 0)])
        S = cs[8:] - cs[:-8]; N = cn[8:] - cn[:-8]                              # window STARTING at grid[k]
        with np.errstate(invalid='ignore', divide='ignore'):
            avg = np.where(N >= 6, S / np.maximum(N, 1), np.nan)               # 3(b) divisor 6/7/8
            sub0 = S / 8.0                                                      # <= 5 h: zeros substituted
        avg = truncate_ppb(np.where(np.isfinite(avg), avg, np.nan))
        sub0 = truncate_ppb(sub0)
        avg = np.where(np.isnan(avg) & (N > 0) & (sub0 > NAAQS_PPB), sub0, avg)  # valid only if > NAAQS
        start = grid[:len(S)]
        for d, day in enumerate(days):
            k = start.get_indexer(day + pd.to_timedelta(np.arange(7, 24), 'h'))   # 3(c) 17 windows
            W = avg[k]
            cnt = np.isfinite(W).sum(0)
            with np.errstate(all='ignore'):
                mx = np.nanmax(np.where(np.isfinite(W), W, -np.inf), 0)
            valid = (cnt >= 13) | ((cnt > 0) & (mx > NAAQS_PPB))                # 3(d)
            out[d, idx[valid]] = mx[valid]
            nval[d, idx] = cnt
    return (days, out, nval) if return_counts else (days, out)


if __name__ == '__main__':
    # self-test on synthetic data
    t = pd.date_range('2022-04-01 01:00', '2022-04-04 00:00', freq='h')        # end-of-hour stamps, UTC
    lst = t - pd.Timedelta(hours=1) - pd.Timedelta(hours=5)                     # hour-beginning EST
    base = np.asarray(40.0 + 20 * np.exp(-((lst.hour - 15) / 3.0) ** 2))        # afternoon peak at 15 LST
    days, m = mda8_epa(base[:, None] + 0.4, t, [-5], '2022-04-01', '2022-04-02')
    print('MDA8 synthetic (EST):', m.ravel(), '(whole ppb)')
    assert np.all(m == np.floor(m)), 'truncation'
    H2 = base.copy(); H2[np.asarray((lst.day == 2) & (lst.hour >= 7) & (lst.hour <= 20))] = np.nan   # long gap on Apr 2
    _, m2, c2 = mda8_epa(H2[:, None], t, [-5], '2022-04-01', '2022-04-02', return_counts=True)
    print('with a 14-h gap on Apr 2: MDA8', m2.ravel(), 'valid windows', c2.ravel())
    assert np.isnan(m2[1, 0]) and c2[1, 0] < 13
    H3 = H2.copy(); H3[np.asarray((lst.day == 2) & (lst.hour == 22))] = 600.0               # one very high hour
    _, m3, c3 = mda8_epa(H3[:, None], t, [-5], '2022-04-01', '2022-04-02', return_counts=True)
    print('gap + high hour: MDA8', m3.ravel(), 'valid windows', c3.ravel(), '(valid via > NAAQS rule)')
    assert m3[1, 0] > 70 and c3[1, 0] < 13
    print('self-test OK')
