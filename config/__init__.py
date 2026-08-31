"""
config — repository configuration package.

Importing this package does two things:

1. re-exports everything from :mod:`config.paths`, so ``from config import ARCHIVE``
   works as well as ``from config.paths import ARCHIVE``;
2. puts the repo's ``functions/`` directory on ``sys.path``, so scripts can
   ``from SE_analysis import get_site_index`` without any path juggling of
   their own.

See :mod:`config.paths` for the full path list and for how to override the
cluster locations with environment variables.
"""

import sys as _sys

from . import paths  # noqa: F401  (re-exported as config.paths)
from .paths import *  # noqa: F401,F403

# Make the repo-local shared functions importable by plain module name.
_functions_dir = str(paths.FUNCTIONS_DIR)
if _functions_dir not in _sys.path:
    _sys.path.insert(0, _functions_dir)
