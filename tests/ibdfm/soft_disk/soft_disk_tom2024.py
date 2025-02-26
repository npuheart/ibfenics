# Copyright (C) 2024 XUAN WANG
#
# This file is part of ibfenics (https://github.com/npuheart/ibfenics)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
#
# email : ibfenics@pengfeima.cn
#
# brief : DISK LID DRIVEN CAVITY

import os
import numpy as np
from loguru import logger
from mshr import *
from fenics import *
from ibfenics1.nssolver import IPCSSolver
from ibfenics1.io import (
    unique_filename,
    create_xdmf_file,
    TimeManager,
    write_paramters,
    write_excel,
)
from local_mesh import *
from ElasticDisk import ElasticDisk
construct_function_space_bc = IPCSSolver.construct_function_space_bc

def output_data(file_fluid, file_solid, u0, p0, f, disp, force, velocity, t, n):
    disp.rename("displacement", "displacement")
    force.rename("force", "force")
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
velocity_d = Function(Vs, name="velocity d")
disp = Function(Vs, name="displacement")
force = Function(Vs, name="force")
disp.interpolate(initial_disp)
disp0 = Function(Vs, name="displacement0")
disp0.interpolate(initial_disp)
ib_interpolation.evaluate_current_points(disp._cpp_object)

# Define fluid solver object
V, Q = construct_function_space_bc(u0, p0)
bcu, bcp = calculate_fluid_boundary_conditions(V, Q)
trial_velocity_solver = IPCSSolver(
    u0, p0, f, dt, nu, stab=stab, alpha=alpha, bcu=bcu, bcp=bcp
)
navier_stokes_solver = IPCSSolver(
    u0, p0, f, dt, nu, stab=stab, alpha=alpha, bcu=bcu, bcp=bcp
)

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
logger.info(f"solid fluid mesh ratio(>2) = {fluid_mesh.hmin() / solid_mesh.hmax()}")
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
    conv=conv,
    delta=delta,
    SAV=SAV,
    nu_s=nu_s,
    n_mesh_fluid=n_mesh_fluid,
    n_mesh_solid=n_mesh_solid,
)

t = dt
time_manager = TimeManager(T, num_steps, 20)
elastic_disk = ElasticDisk(Vs, dt, nu_s, rho)
volume_list = []
for n in range(1, num_steps + 1):
    u_tilde = trial_velocity_solver.solve(bcu)
    # step 1. calculate velocity and pressure
    # logger.info(f"kinematic_energy(u0) {kinematic_energy(u0)}")
    # step 2. interpolate velocity from fluid to solid
    u0_1 = project(u_tilde, Vf_1)
    ib_interpolation.fluid_to_solid(u0_1._cpp_object, velocity._cpp_object)
    # step 3. calculate disp for solid and update current gauss points and dof points

    # step 4. calculate body force.
    elastic_disk.update_velocity(velocity)
    elastic_disk.update_displacement(disp)
    velocity_d = elastic_disk.solve()  # 更新固体的速度
    force.vector()[:] = (velocity_d.vector()[:] - velocity.vector()[:]) / dt 

    # step 5. interpolate force from solid to fluid
    ib_interpolation.solid_to_fluid(f._cpp_object, force._cpp_object)
    u1, p1 = navier_stokes_solver.solve(bcu, bcp)
    navier_stokes_solver.update(u1, p1)
    trial_velocity_solver.update(u1, p1)
    logger.info(f"u0.vector().norm('l2') {u0.vector().norm('l2')}")
    logger.info(f"p0.vector().norm('l2') {p0.vector().norm('l2')}")
    logger.info(f"f.vector().norm('l2') {f.vector().norm('l2')}")
    u0_1 = project(u1, Vf_1)
    ib_interpolation.fluid_to_solid(u0_1._cpp_object, velocity._cpp_object)
    advance_disp_be(disp, velocity, dt)
    elastic_disk.update_displacement(disp)

    ib_interpolation.evaluate_current_points(disp._cpp_object)
    # step 6. update variables and save to file.
    output_data(file_fluid, file_solid, u0, p0, f, disp, force, velocity, t, n)
    volume_list.append(calculate_volume(disp))
    t = n * dt

write_excel(volume_list, file_excel_name)
