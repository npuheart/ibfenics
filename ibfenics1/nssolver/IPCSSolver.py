# Copyright (C) 2024 Pengfei Ma
#
# This file is part of ibfenics1 (https://github.com/npuheart/ibfenics1)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
#
# email : ibfenics1@pengfeima.cn
#
#
# 20241030, a wrapper for the IPCS solver, a solver for the N-S equations

from dolfin import *


class IPCSSolver:
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
        u_ = Function(V)
        # p_n = Function(Q)
        p_ = Function(Q)

        # Define expressions used in variational forms
        U = 0.5 * (u_n + u)
        n = FacetNormal(mesh)
        # f  = Constant((0, 0))
        k = Constant(dt)
        nu = Constant(nu)
        rho = Constant(rho)

        # Define symmetric gradient
        def epsilon(u):
            return sym(nabla_grad(u))

        # Define stress tensor
        def sigma(u, p):
            return 2 * nu * epsilon(u) - p * Identity(len(u))

        # Define variational problem for step 1
        F1 = (
            rho * dot((u - u_n) / k, v) * dx
            + rho * dot(dot(u_n, nabla_grad(u_n)), v) * dx
            + inner(sigma(U, p_n), epsilon(v)) * dx
            + dot(p_n * n, v) * ds
            - dot(nu * nabla_grad(U) * n, v) * ds
            - dot(f, v) * dx
        )
        a1 = lhs(F1)
        L1 = rhs(F1)

        # Define variational problem for step 2
        a2 = dot(nabla_grad(p), nabla_grad(q)) * dx
        L2 = dot(nabla_grad(p_n), nabla_grad(q)) * dx - (1 / k) * div(u_) * q * dx

        # Define variational problem for step 3
        a3 = dot(u, v) * dx
        L3 = dot(u_, v) * dx - k * dot(nabla_grad(p_ - p_n), v) * dx

        # Assemble matrices
        A1 = assemble(a1)
        A2 = assemble(a2)
        A3 = assemble(a3)

        # Apply boundary conditions to matrices
        [bc.apply(A1) for bc in bcu]
        [bc.apply(A2) for bc in bcp]

        self.u_n = u_n
        self.p_n = p_n
        self.u_ = u_
        self.p_ = p_

        self.L1 = L1
        self.L2 = L2
        self.L3 = L3
        self.A1 = A1
        self.A2 = A2
        self.A3 = A3

    def update(self, un, pn):
        self.u_n.assign(un)
        self.p_n.assign(pn)

    def solve(self, bcu, bcp):
        # Step 1: Tentative velocity step
        b1 = assemble(self.L1)
        [bc.apply(b1) for bc in bcu]
        solve(self.A1, self.u_.vector(), b1, "bicgstab", "hypre_amg")

        # Step 2: Pressure correction step
        b2 = assemble(self.L2)
        [bc.apply(b2) for bc in bcp]
        solve(self.A2, self.p_.vector(), b2, "bicgstab", "hypre_amg")

        # Step 3: Velocity correction step
        b3 = assemble(self.L3)
        solve(self.A3, self.u_.vector(), b3, "cg", "sor")
        return self.u_, self.p_
