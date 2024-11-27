# Copyright (C) 2024 Pengfei Ma
#
# This file is part of ibfenics (https://github.com/npuheart/ibfenics)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
#
# email : ibfenics@pengfeima.cn
#
#
# Chorin Solver also called projection method.

from dolfin import *


class ProjectSolver:
    def __init__(self, u0, p0, dt, nu, rho, bcu=None):
        # Reconstruct element space
        mesh = u0.function_space().mesh()
        V = u0.function_space()
        Q = p0.function_space()

        # Define trial and test functions
        u = TrialFunction(V)
        v = TestFunction(V)

        # Define functions for solutions at previous and current time steps
        u1 = Function(V)
        k = Constant(dt)
        nu = Constant(nu)
        rho = Constant(rho)

        # Tentative velocity step
        F1 = (
            rho * (1 / k) * inner(u - u0, v) * dx
            + rho * inner(grad(u0) * u0, v) * dx
            + inner(grad(p0),v)*dx
            + nu * inner(grad(u), grad(v)) * dx
        )
        a1 = lhs(F1)
        L1 = rhs(F1)

        # Assemble matrices
        A1 = assemble(a1)

        # Apply boundary conditions to matrices
        if bcu:
            [bc.apply(A1) for bc in bcu]

        self.u0 = u0
        self.p0 = p0
        self.u1 = u1

        self.L1 = L1
        self.A1 = A1
        self.prec = "amg" if has_krylov_solver_preconditioner("amg") else "default"

    def update(self, u0, p0):
        self.u0.assign(u0)
        self.p0.assign(p0)

    def solve(self, bcu):
        # Compute tentative velocity step
        b1 = assemble(self.L1)
        [bc.apply(self.A1, b1) for bc in bcu]
        solve(self.A1, self.u1.vector(), b1, "bicgstab", "default")
        return self.u1
