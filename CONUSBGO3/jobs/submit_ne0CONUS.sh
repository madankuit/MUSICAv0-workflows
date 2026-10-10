#!/bin/bash
# CONUSBGO3 deliverables for the ne0CONUSne30x8 runs (BASE, noBBCONUS; Apr-Oct 2023 + 2024).
#   w  ne0CONUS -> 1x1 conservative weights (rootxesmf env)
#   c  monitor -> ne0CONUS column indices
#   r  ne30 regression: recompute BASE2022 with the changed regrid code; must equal the file
#   mB, mN  merge BASE / noBBCONUS h2 -> surface O3, one file per year
#   p  point MDA8 + hourly                         (after c, mB, mN)
#   g1 regrid to 1x1  + unified MDA8 (1x1)         (after w, r, mB, mN)
#   g2 regrid to 0.15 + unified MDA8 (0.15)        (after r, mB, mN)
# Usage (from CONUSBGO3/jobs/): bash submit_ne0CONUS.sh <ne30 unified file for the regression> <YYYYMMDD>
set -euo pipefail
REF=${1:?ne30 unified MDA8 file}; CDATE=${2:?YYYYMMDD stamp}
test -f "$REF"; mkdir -p logs
S=bgo3_step.sbatch; E="--experiment ne0CONUS"
w=$(BGO3_ENV=rootxesmf sbatch --parsable -J ne0_weights --time=04:00:00 --mem=64G $S regridding/gen_ne30_to_1x1_weights.py $E)
c=$(sbatch --parsable -J ne0_colidx --time=01:00:00 --mem=16G $S postprocessing/GetMatched_ne30_GivenMonitors_ColumnIndex.py $E)
r=$(sbatch --parsable -J ne30_regress --time=01:00:00 --mem=32G $S regridding/Regrid_ne30_surfO3_to_1x1_conserve.py \
      --scenarios BASE --years 2022 --check-against "$REF" --no-write)
mB=$(sbatch --parsable -J ne0_mergeBASE $S postprocessing/Merge_h2files_hourlysurfO3.py $E --scenarios BASE)
mN=$(sbatch --parsable -J ne0_mergeNOBB $S postprocessing/Merge_h2files_hourlysurfO3.py $E --scenarios noBBCONUS)
p=$(sbatch --parsable -J ne0_points --dependency=afterok:$c:$mB:$mN --mem=64G $S \
      postprocessing/Extract_givenmonitorO3_hourly_dailyMDA8_toNetCDF.py $E)
g1=$(sbatch --parsable -J ne0_regrid1x1 --dependency=afterok:$w:$r:$mB:$mN --mem=64G $S \
      regridding/Regrid_ne30_surfO3_to_1x1_conserve.py $E --target 1x1 --cdate "$CDATE")
g2=$(sbatch --parsable -J ne0_regrid015 --dependency=afterok:$r:$mB:$mN --mem=128G $S \
      regridding/Regrid_ne30_surfO3_to_1x1_conserve.py $E --target 0.15 --cdate "$CDATE")
echo "weights $w | colidx $c | ne30-regress $r | mergeBASE $mB | mergeNOBB $mN | points $p | regrid1x1 $g1 | regrid015 $g2"
