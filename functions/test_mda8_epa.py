#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_mda8_epa.py

Clause-by-clause test of functions/mda8_epa.py against 40 CFR Part 50 Appendix U,
Section 3. `literal()` follows the regulation text line by line (slow, obvious code);
mda8_epa() (vectorized) must agree with it exactly on targeted edge cases and on
random data with random gaps.

  3(a) hours in local standard time; hourly values truncated to 0.001 ppm (whole ppb)
  3(b) 8-h average for every hour, stored at its START hour; >= 6 of 8 hours, divisor
       6 or 7; <= 5 hours: zero substituted, divisor 8, valid only if > NAAQS
       (equivalently: sum of available hours > 0.567 ppm); averages truncated
  3(c) daily max = highest of the 17 averages starting 07:00 ... 23:00 LST
  3(d) daily max valid if >= 13 of 17 averages valid, or if it exceeds the NAAQS

    python test_mda8_epa.py

    18 Oct 2026: VERSION 1.0
"""
import math
import numpy as np, pandas as pd
from mda8_epa import mda8_epa

LEVEL_PPM = 0.070

def literal(hour_begin_lst_ppb, day):
    """hour_begin_lst_ppb: dict {pd.Timestamp hour-beginning LST: ppb float or nan}."""
    def hourly_ppm(t):                                    # 3(a): ppm, 3 decimals, truncated
        v = hour_begin_lst_ppb.get(t, math.nan)
        return math.nan if math.isnan(v) else math.floor(v / 1000.0 * 1000 + 1e-6) / 1000.0
    def eight_hour(start):                                # 3(b)
        vals = [hourly_ppm(start + pd.Timedelta(hours=i)) for i in range(8)]
        have = [v for v in vals if not math.isnan(v)]
        if len(have) >= 6:
            avg = sum(have) / len(have)                   # divisor 6, 7 or 8
        elif sum(have) > 0.567 + 1e-12 and have:          # zero substitution; > NAAQS <=> sum > 0.567 ppm
            avg = sum(have) / 8.0
        else:
            return math.nan
        return math.floor(avg * 1000 + 1e-6) / 1000.0     # truncated to 3 decimals
    avgs = [eight_hour(day + pd.Timedelta(hours=h)) for h in range(7, 24)]     # 3(c): 07..23 starts
    valid = [a for a in avgs if not math.isnan(a)]
    if not valid:
        return math.nan
    mx = max(valid)
    return mx * 1000 if (len(valid) >= 13 or mx > LEVEL_PPM) else math.nan       # 3(d), in ppb

def run_both(lst_hourly, off=-5, days=('2022-07-01', '2022-07-02')):
    """lst_hourly: pd.Series indexed by hour-beginning LST. Build CAM end-of-hour UTC stamps."""
    t_end_utc = lst_hourly.index - pd.Timedelta(hours=off) + pd.Timedelta(hours=1)
    d, m = mda8_epa(lst_hourly.values[:, None], t_end_utc.values, [off], *days)
    ref = [literal(dict(zip(lst_hourly.index, lst_hourly.values)), day) for day in d]
    return m[:, 0], np.array(ref, dtype='float64')

def check(name, got, ref):
    same = np.array_equal(np.isnan(got), np.isnan(ref)) and np.allclose(got[~np.isnan(got)], ref[~np.isnan(ref)], atol=1e-9)
    print(f"{'OK  ' if same else 'FAIL'} {name}: mda8_epa {got.tolist()} | literal {ref.tolist()}")
    assert same, name

idx = pd.date_range('2022-06-30 00:00', '2022-07-04 00:00', freq='h')      # hour-beginning LST
hod = idx.hour.values
base = pd.Series(40 + 30 * np.exp(-((hod - 15) / 3.5) ** 2), index=idx)     # afternoon peak

# 3(a) truncation: 70.999 -> 70, never rounded up
s = pd.Series(np.full(idx.size, 70.999), index=idx); check('3(a) truncation 70.999 -> 70', *run_both(s))
# 3(a) LST and start-hour storage: a single 8-h block of 80 ppb at 07-14 LST only
s = pd.Series(np.full(idx.size, 30.0), index=idx); s[(hod >= 7) & (hod <= 14)] = 80.0
check('3(a)/(c) block 07:00-15:00 LST = first period', *run_both(s))
# same block shifted to 23:00-07:00 next day (last period of the day it starts)
s = pd.Series(np.full(idx.size, 30.0), index=idx); s[(hod >= 23) | (hod <= 6)] = 80.0
check('3(c) 23:00-07:00 block counts for the day it starts', *run_both(s))
# block entirely 00:00-07:00 (not a period start 07..23 of that day -> belongs to previous day's last period)
s = pd.Series(np.full(idx.size, 30.0), index=idx); s[(hod >= 0) & (hod <= 6)] = 80.0
check('3(c) early-morning block only via previous day 23:00 period', *run_both(s))
# 3(b) divisor 6 and 7
for miss in (1, 2):
    s = base.copy(); s[(hod >= 12) & (hod < 12 + miss)] = np.nan
    check(f'3(b) {8-miss} hours -> divisor {8-miss}', *run_both(s))
# 3(b) <= 5 hours, zero substitution, threshold sum > 567 ppb
for v in (113.4, 113.6, 142.0):          # 5 hours: sums 567 / 568 / 710 ppb (after truncation 565/565/710)
    s = pd.Series(np.full(idx.size, np.nan), index=idx); s[(idx.day == 1) & (hod >= 15) & (hod <= 19)] = v
    check(f'3(b) 5 hours of {v} ppb, zero substitution', *run_both(s))
s = pd.Series(np.full(idx.size, np.nan), index=idx); s[(idx.day == 1) & (hod >= 15) & (hod <= 19)] = 113.6
s[(idx.day == 1) & (hod == 15)] = 115.0                                      # sum 115+4*113 = 567 -> not valid
check('3(b) sum exactly 567 ppb -> not valid', *run_both(s))
s[(idx.day == 1) & (hod == 15)] = 116.0                                      # sum 568 -> valid (71 ppb)
check('3(b) sum 568 ppb -> valid, 71 ppb', *run_both(s))
# 3(d) 13 valid periods vs 12. Hours from day1 00:00 + m0 h onward are missing; a period
# starting at s (hours s..s+7) keeps >= 6 hours iff s <= m0 - 6, so starts 07..m0-6 are
# valid: m0 = 25 -> 13 valid (07..19), m0 = 24 -> 12 valid (07..18).
day1 = pd.Timestamp('2022-07-01')
def cut(series, m0):
    s2 = series.copy(); s2[s2.index >= day1 + pd.Timedelta(hours=m0)] = np.nan; return s2
for m0, want in ((25, 'valid'), (24, 'missing')):
    got, ref = run_both(cut(base, m0))
    check(f'3(d) {m0-12} valid periods, max < 70 -> day 1 {want}', got, ref)
    assert (np.isnan(got[0]) == (want == 'missing')), 'case did not exercise the clause'
# 3(d) fewer than 13 valid (12) but max > 70 -> valid; max = 70 after truncation -> missing
for peak, want in ((71.2, 'valid'), (70.9, 'missing')):
    s = pd.Series(np.full(idx.size, 30.0), index=idx); s[(hod >= 9) & (hod <= 16)] = peak
    got, ref = run_both(cut(s, 24))
    check(f'3(d) 12 valid periods, peak {peak} -> day 1 {want}', got, ref)
    assert (np.isnan(got[0]) == (want == 'missing')), 'case did not exercise the clause'
# 3(b) zero substitution must make the day valid only through the 3(d) exceedance rule
s = pd.Series(np.full(idx.size, np.nan), index=idx); s[(idx.day == 1) & (hod >= 15) & (hod <= 19)] = 142.0
got, ref = run_both(s); assert got[0] == 88.0 and ref[0] == 88.0, 'zero-substitution exceedance case'

# random data with random gaps, many cells and offsets
rng = np.random.default_rng(0); nbad = 0; ncmp = 0
for trial in range(200):
    off = int(rng.choice([-5, -6, -7, -8]))
    vals = rng.uniform(5, 120, idx.size)
    vals[rng.random(idx.size) < rng.choice([0.0, 0.1, 0.3, 0.6])] = np.nan
    got, ref = run_both(pd.Series(vals, index=idx), off)
    ncmp += got.size
    nbad += int(not (np.array_equal(np.isnan(got), np.isnan(ref)) and np.allclose(got[~np.isnan(got)], ref[~np.isnan(ref)])))
print(f"{'OK  ' if nbad == 0 else 'FAIL'} random fuzz: {ncmp} daily values, {nbad} trials disagree")
assert nbad == 0
print('ALL APPENDIX U CLAUSE TESTS PASSED')
