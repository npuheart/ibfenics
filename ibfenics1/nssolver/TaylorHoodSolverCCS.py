from dolfin import *
# from mshr import *
import numpy as np


class TaylorHoodSolver:
    def __init__(self, u0, p0, f, dt=0.01, nu=0.01, u_boundary=None):
        # Reconstruct element space
        mesh = u0.function_space().mesh()
        element1 = u0.function_space()._ufl_element
        element2 = p0.function_space()._ufl_element
        TH = element1 * element2
        W = FunctionSpace(mesh, TH)
        N = FacetNormal(mesh)

        # Define variables
        (u, p) = TrialFunctions(W)
        (v, q) = TestFunctions(W)
        k = Constant(dt)
        self.w_ = Function(W)
        self.un, self.pn = Function(W).split(True)

        # Define variational problem
        F = (
            inner((u - self.un) / k, v) * dx
            + inner(grad(self.un) * self.un, v) * dx
            + nu * inner(grad(u), grad(v)) * dx
            - div(v) * p * dx
            - inner(f, v) * dx
        )
        F += q * div(u) * dx

        G = f
        F1 = (
            -inner(grad(p), grad(q)) * dx
            + inner(G, grad(q)) * dx
            - nu * cross(N, self.un) * grad(q) * ds
            - inner(u_boundary, N) * q * ds
        )
        F1 = (
            -inner(grad(p), grad(q)) * dx
            - inner(grad(self.un) * self.un, grad(q)) * dx * dx
        )

        F2 = (
            inner(u, grad(v)) * dx
            + inner(grad(u), grad(v)) * dx
            - 1.0 / nu * inner(G - grad(p), grad(q)) * dx
        )
        F2 = (
            inner(u, grad(v)) * dx
            + inner(grad(u), grad(v)) * dx
            + 1.0 / nu * inner(grad(self.un) * self.un + grad(p), grad(q)) * dx
        )

        a = lhs(F)
        self.W = W
        self.A = assemble(a)
        self.L = rhs(F)

    def update(self, un, pn):
        self.un.assign(un)
        self.pn.assign(pn)

    def solve(self, bcu, bcp):
        # Compute tentative velocity step
        b = assemble(self.L)
        [bc.apply(self.A, b) for bc in bcu]
        [bc.apply(self.A, b) for bc in bcp]
        # solve(self.A, self.w_.vector(), b, "gmres", "hypre_amg")
        solve(self.A, self.w_.vector(), b)

        return self.w_.split(True)


if __name__ == "__main__":
    run_solver()
