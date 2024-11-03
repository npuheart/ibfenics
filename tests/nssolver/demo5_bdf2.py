# Copyright (C) 2024 Pengfei Ma
#
# This file is part of ibfenics (https://github.com/npuheart/ibfenics)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
#
# email : ibfenics@pengfeima.cn
#
# brief : 简单测试流体求解器求解NS方程的正确性(BDF2)


import json
import os
from loguru import logger
import numpy as np
from fenics import *
from mshr import *
from ibfenics.nssolver import TaylorHoodSolverBDF2
from ibfenics.io import unique_filename, create_xdmf_file

note = "demo5"
file_fluid_name = unique_filename(os.path.basename(__file__), note, "/fluid.xdmf")
file_log_name = unique_filename(os.path.basename(__file__), note, "/info.log")
logger.add(file_log_name)


with open("data/demo_navier-stokes_5.json", "r") as file:
    data = json.load(file)

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

    f = Expression((data["f1"], data["f2"]), degree=2, t=0)
    u_exact = Expression((data["u"], data["v"]), degree=2, t=0)
    u_exact_out, p_exact_out = Function(W).split(True)
    p_exact = Expression(data["p"], degree=1, t=0)

    # Define boundary conditions
    bcu = DirichletBC(W.sub(0), u_exact, "on_boundary")
    bcus = [bcu]
    bcps = []

    u0_, _ = Function(W).split(True)
    u0, p0 = Function(W).split(True)
    u_exact.t = 0.0
    u0_.interpolate(u_exact)
    u_exact.t = dt
    u0.interpolate(u_exact)
    navier_stokes_solver = TaylorHoodSolverBDF2(u0_, u0, p0, f, dt, nu)
    for n in range(2, num_steps + 1):
        # 更新时间
        u_exact.t = n * dt
        p_exact.t = n * dt
        f.t = n * dt
        logger.info(f"Step : {n}, Time : {n*dt}, u0(0.5,0.9) : {u0(0.5,0.5)}.")

        navier_stokes_solver.update(u0, p0)
        u1, p1 = navier_stokes_solver.solve(bcus, bcps)
        u0.assign(u1)
        p0.assign(p1)

        # 赋值给u0
        file_fluid.write(u0, n * dt)
        file_fluid.write(p0, n * dt)

        u_exact_out.interpolate(u_exact)
        p_exact_out.interpolate(p_exact)

        print(np.sqrt(assemble(inner(u0 - u_exact, u0 - u_exact) * dx)))
        print(np.sqrt(assemble(inner(p0 - p_exact, p0 - p_exact) * dx)))


if __name__ == "__main__":
    N = 32
    dt = 1.0 / 1024
    T = 0.1
    nu = data["mu"]
    run_solver(dt, nu, T, Nx=N, Ny=N)
