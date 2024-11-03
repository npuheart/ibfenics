# Copyright (C) 2024 Pengfei Ma
#
# This file is part of ibfenics (https://github.com/npuheart/ibfenics)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
#
# email : ibfenics@pengfeima.cn
#
# brief : 简单测试 Taylor-Hood 元求解NS方程的正确性(lid_driven_cavity)

import os
from loguru import logger
import numpy as np
from fenics import *
from mshr import *
from ibfenics.nssolver import TaylorHoodSolver
from ibfenics.io import unique_filename, create_xdmf_file

note = "lid_driven_cavity"
file_fluid_name = unique_filename(os.path.basename(__file__), note, "/fluid.xdmf")
file_log_name = unique_filename(os.path.basename(__file__), note, "/info.log")
logger.add(file_log_name)


# 有了真实解，可以使用零Dirichlet边界条件
def run_solver(dt, nu, T, Nx, Ny):
    num_steps = int(T / dt)
    mesh = UnitSquareMesh(Nx, Ny)
    file_fluid = create_xdmf_file(mesh.mpi_comm(), file_fluid_name)
    logger.info(f"file_fluid_name : {file_fluid_name}, file_log_name : {file_log_name}")
    logger.info(f"mesh size : hmin {mesh.hmin()}, hmax {mesh.hmax()}.")
    logger.info(f"dt = {dt}, nu = {nu}, T = {T}, num_steps = {num_steps}")

    # Define function spaces
    P2 = VectorElement("Lagrange", mesh.ufl_cell(), 2)
    P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
    TH = P2 * P1
    W = FunctionSpace(mesh, TH)

    f = Expression(("0", "0"), degree=2, t=0)
    u_exact = Expression(("1", "0"), degree=2, t=0)
    u_exact_out, p_exact_out = Function(W).split(True)
    p_exact = Expression("0", degree=1, t=0)

    # Define boundary conditions
    upper_flow = Expression(("t", "0"), degree=2, t=0)
    bcu_1 = DirichletBC(W.sub(0), upper_flow, "near(x[1],1.0)")
    bcu_2 = DirichletBC(
        W.sub(0), Constant((0, 0)), "near(x[1],0.0) || near(x[0],0.0) || near(x[0],1.0)"
    )
    bcp = DirichletBC(
        W.sub(1), Constant(0), "near(x[1],0.0) && near(x[0],0.0)", "pointwise"
    )
    bcus = [bcu_1, bcu_2]
    bcps = []

    u0, p0 = Function(W).split(True)
    navier_stokes_solver = TaylorHoodSolver(u0, p0, f, dt, nu)
    for n in range(1, num_steps + 1):
        # 更新时间
        u_exact.t = n * dt
        p_exact.t = n * dt
        upper_flow.t = n * dt
        logger.info(f"Step : {n}, Time : {n*dt}, u0(0.5,0.9) : {u0(0.5,0.5)}.")

        navier_stokes_solver.update(u0, p0)
        u1, p1 = navier_stokes_solver.solve(bcus, bcps)
        u0.assign(u1)
        p0.assign(p1)

        # 赋值给u0
        file_fluid.write(u0, n * dt)
        file_fluid.write(p0, n * dt)

        # TODO: 判断如果存在
        u_exact_out.interpolate(u_exact)
        p_exact_out.interpolate(p_exact)

        print(np.sqrt(assemble(inner(u0 - u_exact, u0 - u_exact) * dx)))
        print(np.sqrt(assemble(inner(p0 - p_exact, p0 - p_exact) * dx)))

        logger.info(f"u0(0.5,0.5) : {u0(0.5,0.5)}.")


if __name__ == "__main__":
    N = 32
    dt = 1.0 / 256
    T = 0.1
    nu = 0.01
    run_solver(dt, nu, T, Nx=N, Ny=N)
