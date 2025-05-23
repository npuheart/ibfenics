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
from scipy.optimize import fsolve


def energy(u, E_hyp, rho):
    return 0.5*rho*assemble(inner(u, u)*dx)  + E_hyp

def modified_energy(u, E_hyp, rho, delta):
    tmp = energy(u, E_hyp, rho) + delta # + 弹性势能
    return tmp

def SAV_X1(un, u1, u2, E_hyp, dt, mu, delta, rho):
    return 2.0/dt*modified_energy(un, E_hyp, rho, delta) + mu*assemble(inner(grad(u2), grad(u2))*dx)

def SAV_X2(un, u1, u2, E_hyp, dt, mu, delta, rho, qn):
    return - 2.0*qn/dt*np.sqrt(modified_energy(un, E_hyp, rho, delta)) + 2.0*mu*assemble(inner(grad(u1), grad(u2))*dx)

def SAV_X3(un, u1, u2, E_hyp, dt, mu, delta, rho):
    return mu*assemble(inner(grad(u1), grad(u1))*dx)

def calculate_SAV(u0, u1, u2, E, dt, nu, delta, qn):
    X1 = SAV_X1(u0, u1, u2, E, dt, nu, delta, 1.0)
    X2 = SAV_X2(u0, u1, u2, E, dt, nu, delta, 1.0, qn)
    X3 = SAV_X3(u0, u1, u2, E, dt, nu, delta, 1.0)
    # 定义目标函数
    def func(x):
        return X1*x*x + X2*x + X3
    # 使用 fsolve 求解
    initial_guess = 1.0  # 初始猜测值
    solution = fsolve(func, initial_guess)
    print(f"解是: {solution[0]}")
    # 直接求解
    # S = (-X2 + np.sqrt(X2*X2-4.0*X1*X3))/(2*X1)
    # if np.isnan(S):
    #     S = -X2/2.0/X1
    return solution[0]

def calculate_SAV_2(u0, u1, u2, E_hyp, dt, nu, delta, qn):
    A0 = delta + assemble(0.5*inner(u1,u1)*dx) + E_hyp
    print(E_hyp)
    A1 = assemble(inner(u1, u2)*dx)
    A2 = 0.5*assemble(inner(u2, u2)*dx)

    B0 = nu*assemble(inner(grad(u1), grad(u1))*dx)
    B1 = 2.0*nu*assemble(inner(grad(u1), grad(u2))*dx)
    B2 = nu*assemble(inner(grad(u2), grad(u2))*dx)
    
    E = lambda x: A0 + A1*x + A2*x*x
    def func(x):
        return 2.0/dt*x*x*E(x) - 2.0*qn/dt*x*np.sqrt(E(x)) + B0 + B1*x + B2*x*x
        # return 2.0/dt*x*x*x*E(x) - 2.0*qn/dt*x*x*np.sqrt(E(x)) + B0*x + B1*x*x + B2*x*x*x
    
    initial_guess = 1.0  # 初始猜测值
    solution = fsolve(func, initial_guess)
    print(f"解是: {solution[0]}")
    return solution[0]


def CAL_SAV_2(
    En = 0.1, delta = 0.1, dt = 0.1, 
    alpha = 0.1, h = 0.1, mu = 0.1, un = 0.1,
    u1 = 0.1 , u2 = 0.1 , qn = 0.1, N = 0.1, 
    w = 0.1, p1 = 0.1, p2 = 0.1, rho = 1.0):
    _1 = 2.0/dt*(En + delta)
    _2 = (alpha*h+mu)*inner(grad(u2),grad(u2))*dx
    _3 = 2.0*qn*np.sqrt(En + delta)/dt
    _4 = (alpha*h+mu)*inner(grad(u1),grad(u2))*dx
    _5 = alpha*h*inner(grad(un),grad(u2))*dx
    _6 = mu*inner(dot(N,grad(u2)), w)*ds
    _7 = inner(N,w)*p2*ds
    _8 = alpha*h*inner(N*grad(u2), w)*ds
    _9 = (alpha*h+mu)*inner(grad(u1),grad(u1))*dx
    _10 = alpha*h*inner(grad(un), grad(u1))*dx 
    _11 = mu*inner(N*grad(u1), w)*ds
    _12 = inner(N,w)*p1*ds
    _13 = alpha*h*inner(N*grad(u1-un), w)*ds
    _14 = 0.5*rho*inner(N*grad(u1)*inner(w,w))*ds


    A = _1 + _2 
    B = -_3 + 2.0*_4 - _5 - _6 + _7 - _8
    C = _9 - _10 - _11 + _12 - _13 - _14

    S = (-B + np.sqrt(B*B - 4.0*A*C))/(2.0*A)
    return S


qn = 1.0
_ = 0

class TaylorHoodSolver_1:
    def __init__(self, u0, p0, dt, nu,stab=False, alpha = 0.1):
        # Reconstruct element space
        mesh = u0.function_space().mesh()
        element1 = u0.function_space()._ufl_element
        element2 = p0.function_space()._ufl_element
        TH = element1*element2
        W = FunctionSpace(mesh, TH)

        # Define variables
        (u, p) = TrialFunctions(W)
        (v, q) = TestFunctions(W)
        k = Constant(dt)
        self.w_ = Function(W)
        self.un, self.pn = Function(W).split(True)

        # Define variational problem
        F = inner((u-self.un)/k, v)*dx + nu * inner(grad(u), grad(v))*dx - div(v)*p * dx
        if stab:
            F += alpha * inner(grad(u-self.un), grad(v))*dx  # 稳定项
        F += q*div(u)*dx
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


class TaylorHoodSolver_2:
    def __init__(self, u0, p0, f, dt, nu, stab=False, alpha = 0.1, conv=True):
        # Reconstruct element space
        mesh = u0.function_space().mesh()
        element1 = u0.function_space()._ufl_element
        element2 = p0.function_space()._ufl_element
        TH = element1*element2
        W = FunctionSpace(mesh, TH)

        # Define variables
        (u, p) = TrialFunctions(W)
        (v, q) = TestFunctions(W)
        k = Constant(dt)
        self.w_ = Function(W)
        self.un, self.pn = Function(W).split(True)

        # Define variational problem
        F = inner(u/k, v)*dx + nu * inner(grad(u), grad(v))*dx + div(v)*p * dx - inner(f, v)*dx
        
        if conv:
            F += inner(grad(self.un)*self.un, v)*dx # 对流项
        
        if stab:
            F += alpha * inner(grad(u), grad(v))*dx  # 稳定项
            
        F += q*div(u)*dx
        a = lhs(F)
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
