#!/usr/bin/env python3
"""
check_paths.py

Verify that the paths declared in config/paths.py actually exist on this
machine, and suggest replacements for the ones that do not.

Unlike its predecessor this script keeps **no path list of its own** — it
introspects config/paths.py, so the config file stays the single source of
truth and the two can never drift apart.

Paths are reported in two groups:

  REPO      files that ship with the repository. These must always exist; a
            miss means the checkout is incomplete (e.g. Git LFS not pulled).
  CLUSTER   data that lives outside the repository. Misses here are normal on
            a machine that does not hold the data, and are what the
            MUSICA_ENV_* environment variables are for.

Usage:
    python check_paths.py                # check everything
    python check_paths.py --cluster      # only the cluster paths
    python check_paths.py --repo         # only the repo-shipped files
    python check_paths.py --suggest      # also hunt for moved files/dirs

MODIFICATION HISTORY:
    VERSION 1.0
    - Initial version
    VERSION 2.0
    - Rewritten to introspect config/paths.py instead of duplicating the list;
      renamed from check_svante_paths.py (the config is no longer Svante-only)
"""

import argparse
import glob
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import paths as P  # noqa: E402


# ============================================================
# COLLECT PATHS FROM THE CONFIG
# ============================================================

def collect():
    """Yield (name, Path, group) for every path constant in config.paths.

    `group` is 'REPO' for anything under the repository, 'CLUSTER' otherwise.
    Dict-valued constants (e.g. TROPOMI_015_DIRS) are expanded per key.
    """
    out = []
    for name in sorted(vars(P)):
        if name.startswith('_') or name.isupper() is False:
            continue
        val = getattr(P, name)
        items = []
        if isinstance(val, Path):
            items = [(name, val)]
        elif isinstance(val, dict) and val and all(isinstance(v, Path) for v in val.values()):
            items = [(f'{name}[{k!r}]', v) for k, v in val.items()]
        for label, p in items:
            group = 'REPO' if _under(p, P.REPO_ROOT) else 'CLUSTER'
            out.append((label, p, group))
    return out


def _under(p, root):
    try:
        p.resolve().relative_to(root.resolve())
        return True
    except ValueError:
        return False


# ============================================================
# SEARCH HELPERS (for --suggest)
# ============================================================

def find_file_candidates(p):
    """Look for a file with the same name, or the same stem with a different
    version tag, in the parent and grandparent directories."""
    candidates = []
    base_stem = re.sub(r'_c\d{8}$', '', p.stem)
    for search_dir in [p.parent, p.parent.parent]:
        if not search_dir.exists():
            continue
        for found in search_dir.rglob(p.name):
            candidates.append((str(found), 'exact name found nearby'))
        for found in glob.glob(str(search_dir / f'*{base_stem}*{p.suffix}')):
            if found not in [c[0] for c in candidates]:
                candidates.append((found, f'similar name (base: {base_stem})'))
    return candidates[:5]


def find_dir_candidates(p):
    """Walk up to the deepest existing ancestor, then look for a directory of
    the same name beneath it; otherwise list the ancestor for context."""
    candidates = []
    ancestor = p
    while not ancestor.exists() and ancestor != ancestor.parent:
        ancestor = ancestor.parent
    if not ancestor.exists() or ancestor == Path('/'):
        return candidates
    for found in ancestor.rglob(p.name):
        if found.is_dir():
            candidates.append((str(found), f"directory '{p.name}' found under {ancestor}"))
    if not candidates and ancestor != p:
        children = [c.name for c in sorted(ancestor.iterdir()) if c.is_dir()][:8]
        if children:
            candidates.append((
                f'Not found. Deepest existing ancestor: {ancestor}\n'
                f'  Contents: {", ".join(children)}',
                'ancestor listing'))
    return candidates[:5]


# ============================================================
# MAIN
# ============================================================

def run(groups, suggest=False):
    entries = [e for e in collect() if e[2] in groups]

    print()
    print('=' * 70)
    print('  MUSICA PATH CHECK')
    print('=' * 70)
    print(f'  config     : {Path(P.__file__).resolve()}')
    print(f'  repo root  : {P.REPO_ROOT}')
    print(f'  data root  : {P.DATA_ROOT}   (MUSICA_ENV_DATA_ROOT)')
    print(f'  home root  : {P.HOME_ROOT}   (MUSICA_ENV_HOME_ROOT)')
    print('=' * 70)

    missing = []
    for group in ('REPO', 'CLUSTER'):
        rows = [e for e in entries if e[2] == group]
        if not rows:
            continue
        print(f'\n  --- {group} ---')
        for name, p, _ in rows:
            if p.exists():
                print(f'  [OK]    {name}')
            else:
                print(f'  [MISS]  {name}')
                print(f'          {p}')
                missing.append((name, p, group))

    n_ok = len(entries) - len(missing)
    print(f'\n{"=" * 70}')
    print(f'  {n_ok} OK  |  {len(missing)} MISSING  (of {len(entries)} checked)')
    print('=' * 70)

    repo_missing = [m for m in missing if m[2] == 'REPO']
    if repo_missing:
        print('\n  WARNING: repo-shipped files are missing. The checkout is')
        print('  incomplete — these should always be present:')
        for name, p, _ in repo_missing:
            print(f'    - {name}: {p}')

    if not missing:
        print('\n  All paths are valid.\n')
        return 0

    if suggest:
        print('\n  Searching for alternatives...\n')
        for name, p, _ in missing:
            print(f'  --- {name} ---')
            print(f'  Missing: {p}')
            cands = find_file_candidates(p) if p.suffix else find_dir_candidates(p)
            if cands:
                print('  Suggestions:')
                for found, reason in cands:
                    print(f'    -> {found}  [{reason}]')
            else:
                print('  No alternatives found automatically.')
            print()
    else:
        print('\n  Re-run with --suggest to hunt for moved files/directories.')
        print('  Cluster misses are expected off-cluster; point the')
        print('  MUSICA_ENV_* variables at your own locations if needed.\n')

    return 1 if repo_missing else 0


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Check the paths declared in config/paths.py.')
    ap.add_argument('--repo', action='store_true', help='only check repo-shipped files')
    ap.add_argument('--cluster', action='store_true', help='only check cluster data paths')
    ap.add_argument('--suggest', action='store_true', help='search for moved files/directories')
    a = ap.parse_args()
    if a.repo and a.cluster:
        groups = {'REPO', 'CLUSTER'}
    elif a.repo:
        groups = {'REPO'}
    elif a.cluster:
        groups = {'CLUSTER'}
    else:
        groups = {'REPO', 'CLUSTER'}
    sys.exit(run(groups, suggest=a.suggest))
