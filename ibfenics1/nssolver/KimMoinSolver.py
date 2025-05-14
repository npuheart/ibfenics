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


class KimMoinSolver:
    def __init__(
        self,
        u_n,
        p_n,
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
        mesh = u_n.function_space().mesh()
        V = u_n.function_space()
        Q = p_n.function_space()

        # Define trial and test functions
        u = TrialFunction(V)
        v = TestFunction(V)
        p = TrialFunction(Q)
        q = TestFunction(Q)

        # Define functions for solutions at previous and current time steps
        # u_n = Function(V)
        u_n_1 = Function(V)
        u_ = Function(V)
        u_hat = 2.0 * u_n - u_n_1
        p_n_1 = Function(Q)
        p_ = Function(Q)

        k = Constant(dt)
        nu = Constant(nu)
        rho = Constant(rho)

        # Tentative velocity step
        F1 = (
            rho*(1/k)*inner(u - u_n, v) * dx + rho * inner(grad(u_hat)*u_hat, v)*dx + \
            nu*inner(0.5*grad(u+u_n), grad(v))*dx - inner(f, v)*dx
        )
        a1 = lhs(F1)
        L1 = rhs(F1)

        # Pressure update
        a2 = inner(grad(p), grad(q)) * dx
        L2 = -(1 / k) * div(u_) * q * dx

        # Velocity update
        a3 = inner(u, v) * dx
        L3 = inner(u_, v) * dx - k * inner(grad(p_), v) * dx

        # Assemble matrices
        A1 = assemble(a1)
        A2 = assemble(a2)
        A3 = assemble(a3)

        # Apply boundary conditions to matrices
        [bc.apply(A1) for bc in bcu]
        [bc.apply(A2) for bc in bcp]

        self.u_n = u_n
        self.u_n_1 = u_n_1
        self.p_n = p_n
        self.u_ = u_
        self.p_ = p_

        self.L1 = L1
        self.L2 = L2
        self.L3 = L3
        self.A1 = A1
        self.A2 = A2
        self.A3 = A3

        self.prec = "amg" if has_krylov_solver_preconditioner("amg") else "default"

    def update(self, u_n, p_n):
        self.u_n_1.assign(self.u_n)  # u_{n-1} = u_n
        self.u_n.assign(u_n)
        self.p_n.assign(p_n)

    def solve(self, bcu, bcp):
        # Compute tentative velocity step
        b1 = assemble(self.L1)
        [bc.apply(self.A1, b1) for bc in bcu]
        solve(self.A1, self.u_.vector(), b1, "bicgstab", "default")

        # Pressure correction
        b2 = assemble(self.L2)
        [bc.apply(self.A2, b2) for bc in bcp]
        [bc.apply(self.p_.vector()) for bc in bcp]
        solve(self.A2, self.p_.vector(), b2, "bicgstab", self.prec)

        # Velocity correction
        b3 = assemble(self.L3)
        [bc.apply(self.A3, b3) for bc in bcu]
        solve(self.A3, self.u_.vector(), b3, "bicgstab", "default")
        # print(assemble(div(self.u_)*dx))
        return self.u_, self.p_
