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
        phi = TrialFunction(Q)
        q = TestFunction(Q)
        # Define functions for solutions at previous and current time steps
        # u0 = Function(V)
        u1 = Function(V)
        u2 = Function(V)
        p1 = Function(Q)
        

        k = Constant(dt)
        nu = Constant(nu)
        rho = Constant(rho)

        # Tentative velocity step
        F1 = (
            rho*(1 / k) * inner(u - u0, v) * dx
            + rho*inner(grad(u0) * u0, v) * dx
            + nu * inner(grad(u), grad(v)) * dx
        )
        a1 = lhs(F1)
        L1 = rhs(F1)

        # project method 
        F2 = rho*inner(u-u1,v)*dx + nu * inner(grad(u), grad(v)) * dx - nu * inner(grad(u1), grad(v)) * dx - inner(f, v) * dx
        a2 = lhs(F2)
        L2 = rhs(F2)
        
        # phi update
        a3 = inner(div(u2),q)*dx
        L3 = inner(div(grad(phi)),q)*dx
         
        # velicity update
        a4 = inner(u2,v)*dx-inner(grad(phi),v)*dx
        L4 = inner(u,v)*dx
        
        #pressure update 
        a5 = p0*q*dx+ rho*phi/k*q*dx-nu*div(grad(phi))*q*dx
        L5 = p*q*dx


        # Assemble matrices
        A1 = assemble(a1)
        A2 = assemble(a2)
        A3 = assemble(a3)
        A4 = assemble(a1)
        A5 = assemble(a2)

        # Apply boundary conditions to matrices
        [bc.apply(A1) for bc in bcu]
        [bc.apply(A5) for bc in bcp]

        self.u0 = u0
        self.p0 = p0
        self.u1 = u1
        self.u2 = u2
        self.p1 = p1

        self.L1 = L1
        self.L2 = L2
        self.L3 = L3
        self.L4 = L4
        self.L5 = L5
        
        
        self.A1 = A1
        self.A2 = A2
        self.A3 = A3
        self.A4 = A4
        self.A5 = A5

        self.prec = "amg" if has_krylov_solver_preconditioner("amg") else "default"

    def update(self, u0, p0):
        self.u0.assign(u0)
        self.p0.assign(p0)

    def solve(self, bcu, bcp):
        # Compute tentative velocity step
        b1 = assemble(self.L1)
        [bc.apply(self.A1, b1) for bc in bcu]
        solve(self.A1, self.u1.vector(), b1, "bicgstab", "default")
        
        b2 = assemble(self.L2)
       
        solve(self.A2, self.phi.vector(), b2, "bicgstab", "default")
        
        b3 = assemble(self.L3)
        solve(self.A3, self.u1.vector(), b3, "bicgstab", "default")
        
        b4 = assemble(self.L4)
        solve(self.A4, self.phi.vector(), b4, "bicgstab", "default")
        
         # Pressure correction
        b5 = assemble(self.L5)
        [bc.apply(self.A5, b5) for bc in bcp]
        [bc.apply(self.p1.vector()) for bc in bcp]
        solve(self.A5, self.p1.vector(), b5, "bicgstab", self.prec)


        
        

        return self.u1, self.p1
    
       