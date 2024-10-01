# Copyright (C) 2024 Pengfei Ma
#
# This file is part of ibfenics (https://github.com/ibfenics)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import numpy
from .cpp import __version__
from .mesh import *
from .interaction import *
from .nssolver import __name__ as nssolver_name

# __all__ = ['__version__', 'mesh', 'interaction', 'nssolver_name']
