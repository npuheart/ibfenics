# Copyright (C) 2024 Pengfei Ma
#
# This file is part of ibfenics1 (https://github.com/npuheart/ibfenics1)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
#
# email : ibfenics1@pengfeima.cn
#
# brief : 简单测试SAV方法求解NS方程的正确性

import os
from loguru import logger
import numpy as np
from fenics import *
from mshr import *
from ibfenics1.nssolver import (
    TaylorHoodSolver_1,
    TaylorHoodSolver_2,
    modified_energy,
    CAL_SAV_2,
)
from ibfenics1.io import unique_filename, create_xdmf_file

note = "lid_driven_cavity"
file_fluid_name = unique_filename(os.path.basename(__file__), note, "/fluid.xdmf")
file_log_name = unique_filename(os.path.basename(__file__), note, "/info.log")
logger.add(file_log_name)


_ = 0  # 弹性势能
delta = 0.1
rho = 1.0
SAV = 1.0
many_sav = []

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
    bcus_1 = [bcu_1, bcu_2]
    bcps_1 = []

    bcu_3 = DirichletBC(W.sub(0), Constant((0, 0)), "on_boundary")
    bcus_2 = [bcu_3]
    bcps_2 = []

    u0, p0 = Function(W).split(True)
    u_old, p_old = Function(W).split(True)
    navier_stokes_solver_1 = TaylorHoodSolver_1(u0, p0, dt, nu)
    navier_stokes_solver_2 = TaylorHoodSolver_2(u0, p0, f, dt, nu)
    for n in range(1, num_steps + 1):
        # 更新时间
        qn = np.sqrt(modified_energy(2 * u_old - u0, _, rho, delta))
        logger.info(f"Step : {n}, Time : {n*dt}, u0(0.5,0.9) : {u0(0.5,0.5)}.")
        logger.info(f"qn : {qn}.")
        upper_flow.t = n * dt
        # 计算两个子问题
        navier_stokes_solver_1.update(u0, p0)
        navier_stokes_solver_2.update(u0, p0)
        u1, p1 = navier_stokes_solver_1.solve(bcus_1, bcps_1)
        u2, p2 = navier_stokes_solver_2.solve(bcus_2, bcps_2)
        # TODO: 构造出最后的解
        # SAV = CAL_SAV_2(2*u_old-u0, u1, u2, _, dt, nu, delta, qn)
        N = FacetNormal(mesh)
        # SAV = CAL_SAV_2(dt, dt, nu, _, delta, _, u0, u1, u2, p0, p1, p2, rho)
        SAV = 1.0
        print(qn, SAV, n)
        u_old.vector()[:] = u0.vector()[:]
        u0.vector()[:] = u1.vector()[:] + SAV * u2.vector()[:]
        p0.vector()[:] = p1.vector()[:] - SAV * p2.vector()[:]
        # 赋值给u0
        file_fluid.write(u0, n * dt)
        file_fluid.write(p0, n * dt)

        print(np.sqrt(assemble(inner(u0, u0) * dx)))
        many_sav.append(SAV)

        logger.info(f"u0(0.5,0.5) : {u0(0.5,0.5)}.")


if __name__ == "__main__":
    N = 32
    dt = 1.0 / 256
    T = 0.1
    nu = 0.01
    run_solver(dt, nu, T, Nx=N, Ny=N)
