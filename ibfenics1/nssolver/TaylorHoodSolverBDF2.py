# Copyright (C) 2024 Pengfei Ma
#
# This file is part of ibfenics1 (https://github.com/npuheart/ibfenics1)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
#
# email : ibfenics1@pengfeima.cn

from dolfin import *
from mshr import *
import numpy as np


class TaylorHoodSolverBDF2:
    @staticmethod
    def construct_function_space_bc(u0, p0):
        W = TaylorHoodSolverBDF2.construct_function_space(u0, p0)
        return W.sub(0), W.sub(1)

    @staticmethod
    def construct_function_space(u0, p0):
        mesh = u0.function_space().mesh()
        element1 = u0.function_space()._ufl_element
        element2 = p0.function_space()._ufl_element
        TH = element1 * element2
        W = FunctionSpace(mesh, TH)
        return W

    def __init__(
        self,
        u0_,
        u0,
        p0,
        f,
        dt,
        nu,
        stab=False,
        alpha=0.1,
        conv=True,
        bdry=None,
        bc_Neumann=None,
    ):
        # Reconstruct element space
        mesh = u0.function_space().mesh()
        W = TaylorHoodSolverBDF2.construct_function_space(u0, p0)

        # Define variables
        (u, p) = TrialFunctions(W)
        (v, q) = TestFunctions(W)
        k = Constant(dt)
        N = FacetNormal(mesh)
        self.w_ = Function(W)
        self.un_, self.un, self.pn = u0_, u0, p0

        # dt approximation
        du_dt = (3.0 * u - 4.0 * self.un + self.un_) / (2.0 * k)

        # Define variational problem
        F = (
            inner(du_dt, v) * dx
            + nu * inner(grad(u), grad(v)) * dx
            - div(v) * p * dx
            - inner(f, v) * dx
        )

        if conv:
            u_hat = 2.0 * self.un - self.un_
            u_grad_u = grad(u_hat) * u_hat
            F += inner(u_grad_u, v) * dx

        if stab:
            F += alpha * inner(grad(3.0 * u - 4.0 * self.un + self.un_), grad(v)) * dx  # 稳定项

        F += q * div(u) * dx
        a = lhs(F)
        self.W = W
        self.A = assemble(a)
        self.L = rhs(F)

    def update(self, un, pn):
        self.un_.assign(self.un)  # u_{n-1} = u_n
        self.un.assign(un)  # u_n = u_{n+1}
        self.pn.assign(pn)  # p_n = p_{n+1}

    def solve(self, bcu, bcp):
        b = assemble(self.L)
        [bc.apply(self.A, b) for bc in bcu]
        [bc.apply(self.A, b) for bc in bcp]
        solve(self.A, self.w_.vector(), b)
        return self.w_.split(True)
