# Copyright (C) 2024 Pengfei Ma
#
# This file is part of ibfenics (https://github.com/npuheart/ibfenics)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
#
# email : ibfenics@pengfeima.cn


from .TaylorHoodSolver import TaylorHoodSolver
from .TaylorHoodSolverBDF2 import TaylorHoodSolverBDF2
from .IPCSSolver import IPCSSolver
from .ChorinSolver import ChorinSolver


class SAVTaylorHoodSolverBDF2:
    from .SAVTaylorHoodSolverBDF2 import TaylorHoodSolverBDF2_1
    from .SAVTaylorHoodSolverBDF2 import TaylorHoodSolverBDF2_2
    from .SAVTaylorHoodSolverBDF2 import modified_energy, calculate_SAV
    from .SAVTaylorHoodSolverBDF2 import (
        construct_function_space,
        construct_function_space_bc,
    )


class SAVTaylorHoodSolver:
    from .SAVTaylorHoodSolver import TaylorHoodSolver_1 as TaylorHoodSolver_1
    from .SAVTaylorHoodSolver import TaylorHoodSolver_2 as TaylorHoodSolver_2
    from .SAVTaylorHoodSolver import modified_energy, CAL_SAV_2
    from .SAVTaylorHoodSolver import (
        construct_function_space,
        construct_function_space_bc,
    )


__name__ = "ibfenics.nssolver"

__all__ = [
    "TaylorHoodSolver",
    "TaylorHoodSolverBDF2",
    "IPCSSolver",
    "ChorinSolver",
    "SAVTaylorHoodSolver",
    "SAVTaylorHoodSolverBDF2",
    "__version__",
    "__name__",
]
