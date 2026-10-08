#!/bin/bash
#SBATCH --job-name=bb_emis_check
#SBATCH --partition=hdr200   # edr nodes have glibc < 2.28 (pyproj ImportError)
#SBATCH --ntasks=1
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%j.log
# Submit from emissions/:  mkdir -p logs && sbatch run_plot_bb_emissions.sbatch [--species CO NO]
set -euo pipefail
source /home/software/anaconda3/2023.07/etc/profile.d/conda.sh && conda activate base
export PYTHONUNBUFFERED=1
python plot_bb_emissions_AprOct_BGO3.py "$@"
