# Copyright (C) 2024 Pengfei Ma
#
# This file is part of ibfenics (https://github.com/npuheart/ibfenics)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
#
# email : ibfenics@pengfeima.cn


from .TaylorHoodSolver import TaylorHoodSolver
from .SAVTaylorHoodSolver import TaylorHoodSolver_1 as TaylorHoodSolver_1
from .SAVTaylorHoodSolver import TaylorHoodSolver_2 as TaylorHoodSolver_2
from .SAVTaylorHoodSolver import modified_energy, CAL_SAV_2
from .TaylorHoodSolverBDF2 import TaylorHoodSolverBDF2

__name__ = "ibfenics.nssolver"

__all__ = [
    'TaylorHoodSolver'    ,
    'TaylorHoodSolver_1'  ,
    'TaylorHoodSolver_2'  ,
    'modified_energy'     ,
    'CAL_SAV_2'           ,
    'TaylorHoodSolverBDF2',
    "__version__"
]


