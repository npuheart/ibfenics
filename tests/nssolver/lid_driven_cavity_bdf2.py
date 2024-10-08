# WARNING: To be modified!
import json
from ibfenics import TaylorHoodSolverBDF2
from fenics import *
from mshr import *
import numpy as np
from datetime import datetime
TaylorHoodSolver = TaylorHoodSolverBDF2.TaylorHoodSolverBDF2

file_id = "data/"+datetime.now().strftime('%Y-%m-%d-%H:%M:%S')


# 有了真实解，可以使用零Dirichlet边界条件
def run_solver(dt, nu, T, Nx, Ny):
    num_steps = int(T / dt)
    mesh = UnitSquareMesh(Nx, Ny)
    print("mesh size : ",mesh.hmin(), mesh.hmax())
    print(f"dt = {dt}, nu = {nu}, T = {T}, num_steps = {num_steps}")

    # Define function spaces
    P2 = VectorElement("Lagrange", mesh.ufl_cell(), 2)
    P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
    TH = P2 * P1
    W = FunctionSpace(mesh, TH)

    f = Expression(("0", "0"), degree=2, t=0)
    u_exact = Expression(("1", "0"), degree=2, t=0)
    # u_exact_out, p_exact_out = Function(W).split(True)
    p_exact = Expression("0", degree=1, t=0)

    # Define boundary conditions
    upper_flow = Expression(("t", "0"), degree=2, t=0)
    bcu_1 = DirichletBC(W.sub(0), upper_flow, "near(x[1],1.0)")
    bcu_2 = DirichletBC(W.sub(0), Constant((0,0)), "near(x[1],0.0) || near(x[0],0.0) || near(x[0],1.0)")
    bcp = DirichletBC(W.sub(1), Constant(0), "near(x[1],0.0) && near(x[0],0.0)", "pointwise")
    bcus = [bcu_1, bcu_2]
    bcps = []

    # output velocity
    ufile = XDMFFile(file_id+"/fluid.xdmf")
    ufile.parameters['rewrite_function_mesh'] = False
    ufile.parameters["functions_share_mesh"] = True
    ufile.parameters["flush_output"] = True

    u0_, p0 = Function(W).split(True)
    u0, p0 = Function(W).split(True)
    # u0.interpolate(u_exact)
    navier_stokes_solver = TaylorHoodSolver(u0_, u0, p0, f, dt, nu)

    for n in range(1, num_steps+1):
        # 更新时间
        # print("Time : ", n*dt)
        print(n, " u0(0.5,0.9) : ", u0(0.5,0.5))
        upper_flow.t = n*dt
        navier_stokes_solver.update(u0, p0)
        u1, p1 = navier_stokes_solver.solve(bcus, bcps)
        u0.assign(u1)
        p0.assign(p1)
        # 赋值给u0
        ufile.write(u0, n*dt)
        ufile.write(p0, n*dt)
        
        print(np.sqrt(assemble(inner(u0,u0)*dx)))


if __name__ == '__main__':
    N = 32
    dt = 0.01
    T = 10.0
    nu = 0.001
    run_solver(dt, nu, T, Nx = N, Ny = N)
