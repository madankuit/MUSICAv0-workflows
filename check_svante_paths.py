"""
check_svante_paths.py

Run this script on Svante to verify all MUSICA data paths are still valid.
For any path that no longer exists, it will search for the most likely
current location and suggest an update.

Usage:
    python check_svante_paths.py
    python check_svante_paths.py --fix      # also writes an updated paths file

MODIFICATION HISTORY:
    VERSION 1.0
    - Initial version
"""

import os
import glob
import argparse
from pathlib import Path

# ============================================================
# ALL PATHS TO CHECK
# Format: (variable_name, path, type)
#   type = 'dir'  → must be an existing directory
#   type = 'file' → must be an existing file
# ============================================================

PATHS = [
    # Root directories
    ("svante_archive",              "/net/fs09/d0/<your_username>/CESM22/archive/",                                   "dir"),
    ("CESM22_root",                 "/net/fs09/d0/<your_username>/CESM22/",                                            "dir"),
    ("Datasets_root",               "/net/fs09/d0/<your_username>/Datasets/",                                          "dir"),
    ("ProcessedData_root",          "/net/fs09/d0/<your_username>/ProcessedData/",                                     "dir"),
    ("ncar_copies_root",            "/net/fs09/d0/<your_username>/ncar_copies/",                                       "dir"),

    # Grid files
    ("SCRIP_ne0CONUSne30x8",        "/net/fs09/d0/<your_username>/CESM22/grids/ne0CONUS_ne30x8_np4_SCRIP.nc",        "file"),
    ("SCRIP_ne30np4",               "/net/fs09/d0/<your_username>/CESM22/grids/ne30np4_091226_pentagons.nc",          "file"),
    ("SCRIP_ne30np4_new",           "/net/fs09/d0/<your_username>/CESM22/grids/ne30np4_grid_c20241105.nc",            "file"),
    ("ESMFmap_015grid_file",        "/net/fs09/d0/<your_username>/CESM22/grids/FV_gridinfo_0.15_c20231204.nc",        "file"),
    ("ESMFmap_01grid_file",         "/net/fs09/d0/<your_username>/CESM22/grids/FV_gridinfo_CAMS_c20210219.nc",        "file"),
    ("ESMFmap_1x1grid_file",        "/net/fs09/d0/<your_username>/CESM22/grids/FV1x1grid_info_c20241105.nc",          "file"),
    ("Regridding_015weights_ne0CONUS", "/net/fs09/d0/<your_username>/CESM22/grids/ne0CONUSne30x8_ESMFmap_0.15x0.15_cubit_conserve_cams_c20231204.nc", "file"),
    ("Regridding_01weights_ne0CONUS",  "/net/fs09/d0/<your_username>/CESM22/grids/ESMFmap_0.1x0.1_ne0CONUSne30x8_cubit_conserve_cams.nc",             "file"),
    ("Regridding_09x125weights_ne30",  "/net/fs09/d0/<your_username>/CESM22/grids/ESMFmap_0.9x1.25_ne30np4_cubit_conserve_cams.nc",                   "file"),
    ("Regridding_1x1weights_ne30",     "/net/fs09/d0/<your_username>/CESM22/grids/ESMFmap_0.9x1.25_ne30np4_cubit_conserve_cams.nc",                   "file"),

    # CESM archive and processed output
    ("processed_output_diri",       "/net/fs09/d0/<your_username>/CESM22/processed_output/",                          "dir"),
    ("July2018_surfh2_merged_diri", "/net/fs09/d0/<your_username>/CESM22/processed_output/July2018_surfh2_merged/",   "dir"),
    ("MUSICATropVCDSave_root",       "/net/fs09/d0/<your_username>/CESM22/Calculated_MUSICA_VCD/TroposphericVCD/",      "dir"),
    ("MUSICAVCDSave_root",          "/net/fs09/d0/<your_username>/CESM22/Calculated_MUSICA_VCD/TotalVCD/",             "dir"),
    ("MetPath",                     "/net/fs09/d0/<your_username>/CESM22/Calculated_MUSICA_VCD/TroposphericVCD/met_f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.1MJuly.TS1base01.nc", "file"),

    # Regridded MUSICA output
    ("Regridded_MUSICA_2018_root",  "/net/fs09/d0/<your_username>/CESM22/Regridded_MUSICA_Output/2018_1330LT_TROPOMIcomp/MassConserve_latlon015_MUSICAoutput/", "dir"),

    # TROPOMI L2
    ("TROPOMI_NO2_L2_diri",         "/net/fs09/d0/<your_username>/Datasets/S5P_L2__NO2____HiR_2/2018/",               "dir"),
    ("TROPOMI_HCHO_L2_diri",        "/net/fs09/d0/<your_username>/Datasets/S5P_L2__HCHO___HiR_2/2018/",               "dir"),
    ("TROPOMI_CO_L2_diri",          "/net/fs09/d0/<your_username>/Datasets/S5P_L2__CO_____HiR_2/2018/",               "dir"),

    # TROPOMI regridded 0.15°
    ("TROPOMI_NO2_015_diri",        "/net/fs09/d0/<your_username>/Datasets/S5P_L2__NO2____HiR_2/regrid_2018_MUSICA015/",  "dir"),
    ("TROPOMI_HCHO_015_diri",       "/net/fs09/d0/<your_username>/Datasets/S5P_L2__HCHO___HiR_2/regrid_2018_MUSICA015/", "dir"),
    ("TROPOMI_CO_015_diri",         "/net/fs09/d0/<your_username>/Datasets/S5P_L2__CO_____HiR_2/regrid_2018_MUSICA015/", "dir"),

    # TROPOMI averaging kernels 0.1°
    ("TROPOMI_NO2AK_01_diri",       "/net/fs09/d0/<your_username>/Datasets/TROPOMI_RegridL2_ToL3/Global/TROPOMI_NO2AK_01deg/",   "dir"),
    ("TROPOMI_HCHOAK_01_diri",      "/net/fs09/d0/<your_username>/Datasets/TROPOMI_RegridL2_ToL3/Global/TROPOMI_HCHOAK_01deg/",  "dir"),
    ("TROPOMI_NO2_01_diri",         "/net/fs09/d0/<your_username>/Datasets/TROPOMI_RegridL2_ToL3/Global/TROPOMI_NO2_01deg/",     "dir"),
    ("TM5SurfaceP_diri",            "/net/fs09/d0/<your_username>/Datasets/TM5MP_Model/",                             "dir"),

    # AQS observations
    ("AQS_diri",                    "/net/fs09/d0/<your_username>/Datasets/AQS/",                                      "dir"),
    ("AQS2018_diri",                "/net/fs09/d0/<your_username>/Datasets/AQS/ForYear2018/",                          "dir"),

    # Emissions
    ("BBEmissions_ne0CONUS_diri",   "/net/fs09/d0/<your_username>/ncar_copies/acom/MUSICA/emissions/qfed2.6_finn/ne0conus30x8/", "dir"),
    ("CAMS_v51_orig_diri",          "/net/fs09/d0/<your_username>/ncar_copies/acom/MUSICA/emissions/cams/CAMS-GLOB-ANTv5.1/CAMS-GLOB-ANT_v5.1_orig/", "dir"),
    ("CAMS_v62_orig_diri",          "/net/fs09/d0/<your_username>/ncar_copies/acom/MUSICA/emissions/cams/CAMS-GLOB-ANT_v6.2/CAMS-GLOB-ANT_v6.2_orig/", "dir"),
    ("NEI2017_01_diri",             "/net/fs09/d0/<your_username>/CESM22/CAMS_withCONUS2017NEI/NEI2017_CONUS_output_01deg/",    "dir"),
    ("NEI2022v2_01_diri",           "/net/fs09/d0/<your_username>/CESM22/CAMS6.2_withCONUS2022v2NEI/NEI2022v2_T1_CONUS_output_01deg/", "dir"),
    ("CAMS_NEI2017_merged_root",    "/net/fs09/d0/<your_username>/CESM22/CAMS_withCONUS2017NEI/GLOB_Merged_CAMS_NEI_01deg/",    "dir"),
    ("CAMS_NEI2022_merged_root",    "/net/fs09/d0/<your_username>/CESM22/CAMS6.2_withCONUS2022v2NEI/globCAMS_conusNEI_01deg/",  "dir"),

    # Project output
    ("DanJaffe_output_diri",        "/net/fs09/d0/<your_username>/ProcessedData/DanJaffeMUSICAPostprocessing/",        "dir"),

    # Local home paths
    ("MUSICA_index_diri",           "/home/<your_username>/Scripts/CESM_analysis/ACP_MUSICANEI_scripts/colidxCSV/",  "dir"),
]


# ============================================================
# SEARCH HELPERS
# ============================================================

def find_file_candidates(path):
    """
    For a missing file path, search in the parent and grandparent directories
    for files with a similar name (same stem, ignoring version suffix).
    Returns a list of (found_path, reason) tuples.
    """
    candidates = []
    p = Path(path)
    stem = p.stem  # e.g. 'FV_gridinfo_0.15_c20231204'
    suffix = p.suffix  # e.g. '.nc'

    # Strip trailing date tag like _c20231204 or _c20241105
    import re
    base_stem = re.sub(r'_c\d{8}$', '', stem)

    for search_dir in [p.parent, p.parent.parent]:
        if not search_dir.exists():
            continue
        # Exact name in different parent
        for found in search_dir.rglob(p.name):
            candidates.append((str(found), "exact name found nearby"))
        # Same base stem, any version date
        pattern = str(search_dir / f"*{base_stem}*{suffix}")
        for found in glob.glob(pattern):
            if found not in [c[0] for c in candidates]:
                candidates.append((found, f"similar name (base: {base_stem})"))

    return candidates[:5]  # cap at 5 suggestions


def find_dir_candidates(path):
    """
    For a missing directory, walk up the path finding the deepest existing
    ancestor, then list its contents to suggest what might have moved.
    Returns a list of (found_path, reason) tuples.
    """
    candidates = []
    p = Path(path)
    dirname = p.name

    # Walk up to find deepest existing ancestor
    ancestor = p
    while not ancestor.exists() and ancestor != ancestor.parent:
        ancestor = ancestor.parent

    if not ancestor.exists() or ancestor == Path('/'):
        return candidates

    # Look for directory with same name under the existing ancestor
    for found in ancestor.rglob(dirname):
        if found.is_dir():
            candidates.append((str(found), f"directory '{dirname}' found under {ancestor}"))

    # If none found, list immediate children of ancestor as context
    if not candidates and ancestor != p:
        children = [str(c) for c in sorted(ancestor.iterdir()) if c.is_dir()][:8]
        if children:
            candidates.append((
                f"Not found. Deepest existing ancestor: {ancestor}\n"
                f"  Contents: {', '.join(Path(c).name for c in children)}",
                "ancestor listing"
            ))

    return candidates[:5]


# ============================================================
# MAIN CHECK
# ============================================================

def run_check(write_fix=False):
    ok = []
    missing = []
    suggestions = {}

    print(f"\n{'='*65}")
    print("  MUSICA PATH CHECK — Svante")
    print(f"{'='*65}\n")

    for varname, path, pathtype in PATHS:
        exists = os.path.isdir(path) if pathtype == 'dir' else os.path.isfile(path)
        if exists:
            ok.append(varname)
            print(f"  [OK]     {varname}")
        else:
            missing.append((varname, path, pathtype))
            print(f"  [MISS]   {varname}")
            print(f"           {path}")

    print(f"\n{'='*65}")
    print(f"  {len(ok)} OK  |  {len(missing)} MISSING")
    print(f"{'='*65}\n")

    if not missing:
        print("  All paths are valid!\n")
        return

    print("  Searching for alternatives...\n")

    for varname, path, pathtype in missing:
        print(f"  --- {varname} ---")
        print(f"  Missing: {path}")

        if pathtype == 'file':
            cands = find_file_candidates(path)
        else:
            cands = find_dir_candidates(path)

        if cands:
            print("  Suggestions:")
            for found, reason in cands:
                print(f"    → {found}  [{reason}]")
            suggestions[varname] = cands[0][0]  # best guess = first candidate
        else:
            print("  No alternatives found automatically.")
            suggestions[varname] = None
        print()

    # Optionally write updated paths file
    if write_fix:
        _write_updated_paths(suggestions, missing)


def _write_updated_paths(suggestions, missing):
    """
    Append an UPDATE section to svante_MUSICA_paths.py with suggested replacements.
    """
    outfile = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'svante_MUSICA_paths_updated.py')
    missing_varnames = {v for v, _, _ in missing}

    # Re-read original file
    orig = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'svante_MUSICA_paths.py')
    with open(orig) as f:
        lines = f.readlines()

    updated_lines = []
    for line in lines:
        replaced = False
        for varname, path, pathtype in missing:
            if line.strip().startswith(varname + ' ') and '=' in line and suggestions.get(varname):
                new_path = suggestions[varname]
                indent = len(line) - len(line.lstrip())
                updated_lines.append(f"{' '*indent}{varname:<40} = '{new_path}'  # updated by check_svante_paths.py\n")
                replaced = True
                break
        if not replaced:
            updated_lines.append(line)

    with open(outfile, 'w') as f:
        f.writelines(updated_lines)

    print(f"  Updated paths written to: {outfile}")
    print("  Review the suggestions and rename to svante_MUSICA_paths.py when satisfied.\n")


# ============================================================
# ENTRY POINT
# ============================================================

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Check MUSICA data paths on Svante.')
    parser.add_argument('--fix', action='store_true',
                        help='Write svante_MUSICA_paths_updated.py with suggested replacements')
    args = parser.parse_args()
    run_check(write_fix=args.fix)
