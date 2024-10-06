from dolfin import *
from mshr import *
import numpy as np

class TaylorHoodSolverBDF2:
    def __init__(self, u0_, u0, p0, f, dt, nu, stab=False, alpha = 0.1):
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
        self.un , self.pn = Function(W).split(True)
        self.un_, self.pn_ = Function(W).split(True)
        
        # dt approximation
        du_dt = (3.0*u-4.0*self.un+self.un_)/(2.0*k)
        u_hat = 2.0*self.un - self.un_
        u_grad_u = grad(u_hat)*u_hat

        # Define variational problem
        F = inner(du_dt, v)*dx + inner(u_grad_u, v)*dx + \
            nu * inner(grad(u), grad(v))*dx - div(v)*p * dx - inner(f, v)*dx
        if stab:
            F += alpha * inner(grad(u-self.un), grad(v))*dx  # 稳定项
        F += q*div(u)*dx
        a = lhs(F)
        self.W = W
        self.A = assemble(a)
        self.L = rhs(F)
        self.un.assign(u0_)

    def update(self, un, pn):
        self.un_.assign(self.un) # u_{n-1} = u_n
        self.un.assign(un)       # u_n = u_{n+1}
        self.pn.assign(pn)       # p_n = p_{n+1}

    def solve(self, bcu, bcp):
        b = assemble(self.L)
        [bc.apply(self.A, b) for bc in bcu]
        [bc.apply(self.A, b) for bc in bcp]
        solve(self.A, self.w_.vector(), b)
        return self.w_.split(True)

