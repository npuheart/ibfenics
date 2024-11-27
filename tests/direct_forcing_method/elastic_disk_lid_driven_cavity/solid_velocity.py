
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

class SolidVelocity:
    def __init__(self, w0, u0, p0, dt, rho_s, nu_s, bc=None):
        self.cell = triangle
        P2 = VectorElement("Lagrange", self.cell, 2)
        P1 = FiniteElement("Lagrange", self.cell, 1)
        R = FiniteElement("Real", self.cell, 0)
        ME = FunctionSpace(mesh, MixedElement([P2, P2, P1]))

        test_state = TestFunction(ME)
        dstate     = TrialFunction(ME)  

        # Define functions for solutions at previous and current time steps
        state      = Function(ME)  # current solution
        state0     = Function(ME)  # solution from previous converged step

        k = Constant(dt)
        rho_s = Constant(rho_s)
        nu_s = Constant(nu_s)
        f = Function(ME.sub(0).collapse())
        
        x,  v,  q  = split(test_state)
        dw, du, dp = split(dstate)
        w,  u,  p  = split(state)
        w0, u0, p0 = split(state0)
        
        def calculate_NH_model():  # 这段加入中
            F = grad(state.sub(1))
            H = det(F)*inv(F).T
            P = nu_s * (F - inv(F).T) - p*H
            return P
        
        F = grad(state.sub(1))
        H = det(F)*inv(F).T
        F0 = inner(w-w0,x)*dx - k*inner(u,x)*dx
        F1 = rho_s*inner(u-u0,v)*dx + k*inner(calculate_NH_model(), grad(v))*dx - rho_s*inner(f,v)*dx
        F2 = k*inner(H,grad(u))*q*dx
        F  = F0 + F1 + F2
        
        if bc:
            [bc.apply(F) for bc in bc]

        self.state = state
        self.J = derivative(F, self.state, dstate)
        self.F = F
        self.ME = ME
        
    def update(self, state0):
        self.state0.assign(state0)

    def solve(self):
        solve(self.F == 0, self.state, J=self.J, solver_parameters={"newton_solver": {
        "linear_solver": "mumps", "absolute_tolerance": 1.0e-7, "relative_tolerance": 1e-7}})
        return self.state

if __name__ == "__main__":
    mesh = UnitSquareMesh(8, 8)

    # 创建初始条件和边界条件
    w0 = Expression(("0.0", "0.0"), degree=1)
    u0 = Expression(("0.0", "0.0"), degree=1)   
    p0 = Constant(0)
    dt = 0.01
    rho_s = Constant(1.0)
    nu_s = Constant(0.5)

# 创建边界条件
    def boundary(x, on_boundary):
        return on_boundary
    
    ME = SolidVelocity(w0, u0, p0, dt, rho_s, nu_s).ME 

    bc = DirichletBC(ME.sub(0).collapse(), u0, boundary)

    sv = SolidVelocity(w0, u0, p0, dt, rho_s, nu_s)

    # 时间步进
    t = 0
    T = 1.0
    while t < T:
        
        state = sv.solve()
        t += dt
        SolidVelocity.update()

    # 绘制最终的解
    plot(state.sub(0).collapse())




