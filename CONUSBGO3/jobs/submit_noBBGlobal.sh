#!/bin/bash
# Add the noBBGlobal scenario to the CONUSBGO3 deliverables (Oct 2026).
#   1 merge   h2 -> surface O3, one file per year           (noBBGlobal only)
#   2 points  monitor MDA8 + hourly files                   (after 1)
#   3 check   recompute BASE2022 with the new code; must equal the delivered file
#   4 regrid  noBBGlobal hourly 1x1 + unified 4-scenario MDA8, the other three
#             scenarios copied from the delivered file      (after 1 and 3)
# Usage (from CONUSBGO3/jobs/): bash submit_noBBGlobal.sh <delivered unified file> <YYYYMMDD>
set -euo pipefail
OLD=${1:?delivered 3-scenario unified MDA8 file}; CDATE=${2:?YYYYMMDD stamp}
test -f "$OLD"
mkdir -p logs
j1=$(sbatch --parsable -J bgo3_merge  bgo3_step.sbatch postprocessing/Merge_h2files_hourlysurfO3.py --scenarios noBBGlobal)
j2=$(sbatch --parsable -J bgo3_points --dependency=afterok:$j1 bgo3_step.sbatch \
       postprocessing/Extract_givenmonitorO3_hourly_dailyMDA8_toNetCDF.py --scenarios noBBGlobal)
j3=$(sbatch --parsable -J bgo3_check  bgo3_step.sbatch regridding/Regrid_ne30_surfO3_to_1x1_conserve.py \
       --scenarios BASE --years 2022 --check-against "$OLD" --no-write)
j4=$(sbatch --parsable -J bgo3_regrid --dependency=afterok:$j1:$j3 bgo3_step.sbatch \
       regridding/Regrid_ne30_surfO3_to_1x1_conserve.py --scenarios noBBGlobal --append-to "$OLD" --cdate "$CDATE")
echo "merge $j1 | points $j2 | check $j3 | regrid $j4"
