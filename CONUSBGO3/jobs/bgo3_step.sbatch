#!/bin/bash
#SBATCH --partition=hdr200        # edr nodes have glibc < 2.28 (pyproj/cartopy ImportError)
#SBATCH --ntasks=1
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%j.log
# One CONUSBGO3 pipeline step: sbatch -J <name> bgo3_step.sbatch <script relative to CONUSBGO3/> [args]
# Submitted by submit_noBBGlobal.sh; run from CONUSBGO3/jobs/.
set -euo pipefail
source /home/software/anaconda3/2023.07/etc/profile.d/conda.sh && conda activate base
export PYTHONUNBUFFERED=1
cd "$SLURM_SUBMIT_DIR/.."
echo "host $(hostname)  git $(git rev-parse --short HEAD)  $(date)"
python "$@"
