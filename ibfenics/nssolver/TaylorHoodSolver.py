# Copyright (C) 2024 Pengfei Ma
#
# This file is part of ibfenics (https://github.com/npuheart/ibfenics)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
#
# email : ibfenics@pengfeima.cn

from dolfin import *
from mshr import *
import numpy as np


class TaylorHoodSolver:
    def __init__(
        self,
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
        element1 = u0.function_space()._ufl_element
        element2 = p0.function_space()._ufl_element
        TH = element1 * element2
        W = FunctionSpace(mesh, TH)

        # Define variables
        (u, p) = TrialFunctions(W)
        (v, q) = TestFunctions(W)
        k = Constant(dt)
        N = FacetNormal(mesh)
        self.w_ = Function(W)
        self.un, self.pn = Function(W).split(True)

        if bdry is not None:
            ds = Measure("ds", domain=mesh, subdomain_data=bdry)

        # Define variational problem
        F = (
            inner((u - self.un) / k, v) * dx
            + nu * inner(grad(u), grad(v)) * dx
            - div(v) * p * dx
            - inner(f, v) * dx
        )

        if bdry is not None:
            for i in bc_Neumann:
                F -= inner(nu * dot(N, grad(u)) - p * N, v) * ds(i)

        if conv:
            F += inner(grad(self.un) * self.un, v) * dx

        if stab:
            F += alpha * inner(grad(u - self.un), grad(v)) * dx  # 稳定项

        F += q * div(u) * dx
        a = lhs(F)
        self.W = W
        self.A = assemble(a)
        self.L = rhs(F)

    def update(self, un, pn):
        self.un.assign(un)
        self.pn.assign(pn)

    def solve(self, bcu, bcp):
        b = assemble(self.L)
        [bc.apply(self.A, b) for bc in bcu]
        [bc.apply(self.A, b) for bc in bcp]
        solve(self.A, self.w_.vector(), b)
        return self.w_.split(True)
