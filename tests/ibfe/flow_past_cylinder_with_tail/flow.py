# Copyright (C) 2024 Pengfei Ma
#
# This file is part of ibfenics1 (https://github.com/npuheart/ibfenics1)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
#
# email : ibfenics1@pengfeima.cn
#
# brief : 简单测试 Taylor-Hood 元求解NS方程的正确性(圆柱绕流问题)

import os
from loguru import logger
import numpy as np
from fenics import *
# from mshr import *
from ibfenics1.nssolver import TaylorHoodSolver
from ibfenics1.io import unique_filename, create_xdmf_file

note = "lid_driven_cavity"
file_fluid_name = unique_filename(os.path.basename(__file__), note, "/fluid.xdmf")
file_log_name = unique_filename(os.path.basename(__file__), note, "/info.log")
logger.add(file_log_name)


# 有了真实解，可以使用零Dirichlet边界条件
def run_solver(dt, nu, T, Nx, Ny):
    # Create mesh
    channel = Rectangle(Point(0, 0), Point(2.2, 0.41))
    cylinder = Circle(Point(0.2, 0.2), 0.05)
    domain = channel - cylinder
    mesh = generate_mesh(domain, 64)

    # Set parameters
    num_steps = int(T / dt)
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

    # Define boundary conditions
    inflow = "near(x[0], 0)"
    outflow = "near(x[0], 2.2)"
    walls = "near(x[1], 0) || near(x[1], 0.41)"
    cylinder = "on_boundary && x[0]>0.1 && x[0]<0.3 && x[1]>0.1 && x[1]<0.3"
    inflow_profile = ("4.0*1.5*x[1]*(0.41 - x[1]) / pow(0.41, 2)", "0")
    bcu_inflow = DirichletBC(W.sub(0), Expression(inflow_profile, degree=2), inflow)
    bcu_walls = DirichletBC(W.sub(0), Constant((0, 0)), walls)
    bcu_cylinder = DirichletBC(W.sub(0), Constant((0, 0)), cylinder)
    bcp_outflow = DirichletBC(W.sub(1), Constant(0), outflow)
    bcus = [bcu_inflow, bcu_walls, bcu_cylinder]
    bcps = [bcp_outflow]

    u0, p0 = Function(W).split(True)
    navier_stokes_solver = TaylorHoodSolver(u0, p0, f, dt, nu)
    for n in range(1, num_steps + 1):
        # 更新时间
        # upper_flow.t = n*dt
        logger.info(f"Step : {n}, Time : {n*dt}, u0(0.1,0.1) : {u0(0.1,0.1)}.")

        navier_stokes_solver.update(u0, p0)
        u1, p1 = navier_stokes_solver.solve(bcus, bcps)
        u0.assign(u1)
        p0.assign(p1)

        # 赋值给u0
        file_fluid.write(u0, n * dt)
        file_fluid.write(p0, n * dt)


if __name__ == "__main__":
    N = 32
    dt = 1.0 / 1000.0
    T = 10
    nu = 0.001
    run_solver(dt, nu, T, Nx=N, Ny=N)
