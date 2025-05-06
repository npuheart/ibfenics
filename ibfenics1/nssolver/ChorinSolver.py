# Copyright (C) 2024 Pengfei Ma
#
# This file is part of ibfenics1 (https://github.com/npuheart/ibfenics1)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
#
# email : mapengfei@mail.nwpu.edu.cn
#
#
# Chorin Solver also called projection method.

from dolfin import *


class ChorinSolver:
    def __init__(
        self,
        u0,
        p0,
        f,
        dt,
        nu,
        rho=1.0,
        stab=False,
        alpha=0.1,
        conv=True,
        bdry=None,
        bc_Neumann=None,
        bcu=None,
        bcp=None,
    ):
        # Reconstruct element space
        mesh = u0.function_space().mesh()
        V = u0.function_space()
        Q = p0.function_space()

        # Define trial and test functions
        u = TrialFunction(V)
        v = TestFunction(V)
        p = TrialFunction(Q)
        q = TestFunction(Q)

        # Define functions for solutions at previous and current time steps
        # u0 = Function(V)
        u1 = Function(V)
        p1 = Function(Q)

        k = Constant(dt)
        nu = Constant(nu)
        rho = Constant(rho)

        # Tentative velocity step
        F1 = (
            (1 / k) * inner(u - u0, v) * dx
            + nu * inner(grad(u), grad(v)) * dx
            - inner(f, v) * dx
        )
        a1 = lhs(F1)
        L1 = rhs(F1)
        
        if conv:
            F1 += inner(grad(u0) * u0, v) * dx
            

        # Pressure update
        a2 = inner(grad(p), grad(q)) * dx
        L2 = -(1 / k) * div(u1) * q * dx

        # Velocity update
        a3 = inner(u, v) * dx
        L3 = inner(u1, v) * dx - k * inner(grad(p1), v) * dx

        # Assemble matrices
        A1 = assemble(a1)
        A2 = assemble(a2)
        A3 = assemble(a3)

        # Apply boundary conditions to matrices
        [bc.apply(A1) for bc in bcu]
        [bc.apply(A2) for bc in bcp]

        self.u0 = u0
        self.p0 = p0
        self.u1 = u1
        self.p1 = p1

        self.L1 = L1
        self.L2 = L2
        self.L3 = L3
        self.A1 = A1
        self.A2 = A2
        self.A3 = A3

        self.prec = "amg" if has_krylov_solver_preconditioner("amg") else "default"

    def update(self, u0, p0):
        self.u0.assign(u0)
        self.p0.assign(p0)

    def solve(self, bcu, bcp):
        # Compute tentative velocity step
        b1 = assemble(self.L1)
        [bc.apply(self.A1, b1) for bc in bcu]
        solve(self.A1, self.u1.vector(), b1, "bicgstab", "default")

        # Pressure correction
        b2 = assemble(self.L2)
        [bc.apply(self.A2, b2) for bc in bcp]
        [bc.apply(self.p1.vector()) for bc in bcp]
        solve(self.A2, self.p1.vector(), b2, "bicgstab", self.prec)

        # Velocity correction
        b3 = assemble(self.L3)
        [bc.apply(self.A3, b3) for bc in bcu]
        solve(self.A3, self.u1.vector(), b3, "bicgstab", "default")
        # print(assemble(div(self.u1)*dx))
        return self.u1, self.p1
