# Copyright (C) 2024 Pengfei Ma
#
# This file is part of ibfenics1 (https://github.com/npuheart/ibfenics1)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
#
# email : mapengfei@mail.nwpu.edu.cn
#
# brief : 测试IBFE方法的正确性

import os
import numpy as np
from loguru import logger
# from mshr import *
from fenics import *
from ibfenics1.nssolver import SAVTaylorHoodSolver
from ibfenics1.io import (
    unique_filename,
    create_xdmf_file,
    TimeManager,
    write_paramters,
    write_excel,
    write_excel_sheets,
)
from local_mesh import *

TaylorHoodSolver_1 = SAVTaylorHoodSolver.TaylorHoodSolver_1
TaylorHoodSolver_2 = SAVTaylorHoodSolver.TaylorHoodSolver_2
modified_energy = SAVTaylorHoodSolver.modified_energy
CAL_SAV = SAVTaylorHoodSolver.CAL_SAV_2
construct_function_space_bc = SAVTaylorHoodSolver.construct_function_space_bc

# Define boundary conditions for fluid solver
def calculate_fluid_boundary_conditions(V, Q):
    bcu_1 = DirichletBC(V, Constant((0, 0)), "near(x[1],1.0)")
    bcu_2 = DirichletBC(
        V, Constant((0, 0)), "near(x[1],0.0) || near(x[0],0.0) || near(x[0],1.0)"
    )
    bcp_1 = DirichletBC(Q, Constant(0), "near(x[1],0.0) && near(x[0],0.0)", "pointwise")
    bcu = [bcu_1, bcu_2]
    bcp = [bcp_1]
    return bcu, bcp


def calculate_fluid_boundary_conditions_sav(V, Q):
    bcu_1 = DirichletBC(V, Constant((0, 0)), "on_boundary")
    bcp_1 = DirichletBC(Q, Constant(0), "near(x[1],0.0) && near(x[0],0.0)", "pointwise")
    bcu = [bcu_1]
    bcp = [bcp_1]
    return bcu, bcp


# Define solid constituitive model
def calculate_constituitive_model(disp, vs, us):
    F = grad(disp)
    P = nu_s * (F - inv(F).T)
    F2 = inner(P, grad(vs)) * dx + inner(us, vs) * dx
    a2 = lhs(F2)
    L2 = rhs(F2)
    A2 = assemble(a2)
    return A2, L2


def output_data(file_fluid, file_solid, u0, p0, f, disp, force, velocity, t, n):
    if time_manager.should_output(n):
        logger.info(f"time: {t}, step: {n}, output...")
        file_fluid.write(u0, t)
        file_fluid.write(p0, t)
        file_fluid.write(f, t)
        file_solid.write(disp, t)
        file_solid.write(force, t)
        file_solid.write(velocity, t)


# Create functions for fluid
u0 = Function(Vf, name="velocity")
u0_1 = Function(Vf_1, name="velocity 1st order")
p0 = Function(Vp, name="pressure")
f = Function(Vf_1, name="force")

# Create functions for solid
velocity = Function(Vs, name="velocity")
disp = Function(Vs, name="displacement")
force = Function(Vs, name="force")
disp.interpolate(initial_disp)
ib_interpolation.evaluate_current_points(disp._cpp_object)

# Define fluid solver object
V, Q = construct_function_space_bc(u0, p0)
bcus_1, bcps_1 = calculate_fluid_boundary_conditions(V, Q)
bcus_2, bcps_2 = calculate_fluid_boundary_conditions_sav(V, Q)
navier_stokes_solver_1 = TaylorHoodSolver_1(u0, p0, dt, nu, stab=stab, alpha=alpha)
navier_stokes_solver_2 = TaylorHoodSolver_2(
    u0, p0, f, dt, nu, stab=stab, alpha=alpha, conv=conv
)

# Define trial and test functions for solid solver
us = TrialFunction(Vs)
vs = TestFunction(Vs)
A2, L2 = calculate_constituitive_model(disp, vs, us)

# Define output path
file_log_name = unique_filename(os.path.basename(__file__), str(dt), "/info.log")
file_solid_name = unique_filename(os.path.basename(__file__), str(dt), "/solid.xdmf")
file_fluid_name = unique_filename(os.path.basename(__file__), str(dt), "/fluid.xdmf")
file_excel_name = unique_filename(os.path.basename(__file__), str(dt), "/volume.xlsx")
file_param_name = unique_filename(
    os.path.basename(__file__), str(dt), "/parameters.json"
)
file_solid = create_xdmf_file(solid_mesh.mpi_comm(), file_solid_name)
file_fluid = create_xdmf_file(fluid_mesh.mpi_comm(), file_fluid_name)
logger.add(file_log_name)
logger.info(f"file_solid_name: {file_solid_name}, file_fluid_name: {file_fluid_name}")
logger.info(f"file_excel_name: {file_excel_name}, file_log_name: {file_log_name}")
logger.info(f"solid_mesh.hmax() {solid_mesh.hmax()}, hmin() {solid_mesh.hmin()}")
logger.info(f"fluid_mesh.hmax() {fluid_mesh.hmax()}, hmin() {fluid_mesh.hmin()}")
logger.info("solid fluid mesh ratio(>2) = ", fluid_mesh.hmin() / solid_mesh.hmax())
logger.info(
    f"fluid_mesh.num_cells() {fluid_mesh.num_cells()}, Vp.dim() {Vp.dim()}, Vf.dim() {Vf.dim()}"
)
logger.info(f"solid_mesh.num_cells() {solid_mesh.num_cells()}, Vs.dim() {Vs.dim()}")
write_paramters(
    file_param_name,
    T=T,
    dt=dt,
    num_steps=num_steps,
    rho=rho,
    nu=nu,
    alpha=alpha,
    stab=stab,
    delta=delta,
    SAV=SAV,
    nu_s=nu_s,
    n_mesh_fluid=n_mesh_fluid,
    n_mesh_solid=n_mesh_solid,
)

t = dt
time_manager = TimeManager(T, num_steps, 20)
qn = np.sqrt(total_energy(u0, disp)+delta)
volume_list = []
for n in range(1, num_steps + 1):
    # step 1. calculate velocity and pressure
    En = total_energy(u0, disp)
    u1, p1 = navier_stokes_solver_1.solve(bcus_1, bcps_1)
    u2, p2 = navier_stokes_solver_2.solve(bcus_2, bcps_2)
    SAV = CAL_SAV(En, delta, dt, alpha, h, nu, u0, u1, u2, qn, N, Function(Vf), p1, p2, rho)
    S = CAL_SAV(total_energy(u0, disp), delta, dt, alpha, h, nu, u0, u1, u2, qn, FacetNormal(fluid_mesh), Function(Vf), p1, p2, rho)
    R = np.sqrt(total_energy(u0, disp)+delta)
    Q = S*R
    qn = Q
    SAV = S
    logger.info(f"S                      {S}")
    logger.info(f"R                      {R}")
    logger.info(f"Q                      {Q}")
    u0.vector()[:] = u1.vector()[:] + SAV * u2.vector()[:]
    p0.vector()[:] = p1.vector()[:] - SAV * p2.vector()[:]
    navier_stokes_solver_1.update(u0, p0)
    navier_stokes_solver_2.update(u0, p0)
    logger.info(f"u0.vector().norm('l2') {u0.vector().norm('l2')}")
    logger.info(f"p0.vector().norm('l2') {p0.vector().norm('l2')}")
    logger.info(f"f.vector().norm('l2') {f.vector().norm('l2')}")
    logger.info(f"kinematic_energy(u0) {kinematic_energy(u0)}")
    logger.info(f"total_energy(u0)       {total_energy(u0, disp)}")
    # step 2. interpolate velocity from fluid to solid
    u0_1 = project(u0, Vf_1)
    ib_interpolation.fluid_to_solid(u0_1._cpp_object, velocity._cpp_object)
    # step 3. calculate disp for solid and update current gauss points and dof points
    advance_disp_be(disp, velocity, dt)
    ib_interpolation.evaluate_current_points(disp._cpp_object)
    # step 4. calculate body force.
    b2 = assemble(L2)
    solve(A2, force.vector(), b2)
    # step 5. interpolate force from solid to fluid
    ib_interpolation.solid_to_fluid(f._cpp_object, force._cpp_object)
    # step 6. update variables and save to file.
    output_data(file_fluid, file_solid, u0, p0, f, disp, force, velocity, t, n)
    volume_list.append(calculate_volume(disp))
    t = n * dt

write_excel(volume_list, file_excel_name)
