# Copyright (C) 2024 Pengfei Ma
#
# This file is part of ibfenics (https://github.com/npuheart/ibfenics)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
#
# email : ibfenics@pengfeima.cn
#
# brief : 测试IBFE方法的正确性

import os
import numpy as np
from loguru import logger

from mshr import *
from fenics import *

from ibfenics.nssolver import TaylorHoodSolver
from ibfenics.io import unique_filename, create_xdmf_file, TimeManager, write_paramters

from ref_coordinates import FiberForce
from local_mesh import *
from post_processing import *

# Define time parameters
T =  0.0000005
dt = 0.0000001
num_steps = int(T/dt)
time_manager = TimeManager(T, num_steps, 20)

# Define fluid parameters
nu = 0.01

# Define stablization parameters
alpha = 1.0*dt
stab  = False
delta = 0.1
SAV   = 1.0

def advance_disp_be(disp, velocity, dt):
    disp.vector()[:] = velocity.vector()[:]*dt + disp.vector()[:]

# Define boundary conditions for fluid solver
def calculate_fluid_boundary_conditions(W):
    bcu_1 = DirichletBC(W.sub(0), Constant((0,0)), "near(x[1],1.0)")
    bcu_2 = DirichletBC(W.sub(0), Constant((0,0)), "near(x[1],0.0) || near(x[0],0.0) || near(x[0],1.0)")
    bcp_1 = DirichletBC(W.sub(1), Constant(0), "near(x[1],0.0) && near(x[0],0.0)", "pointwise")
    bcu = [bcu_1, bcu_2]
    bcp = [bcp_1]
    return bcu, bcp

# TODO: Define solid constituitive model
def calculate_constituitive_model(disp, vs, us):
    fiber_force = FiberForce()
    F2 = -inner(fiber_force, vs)*dx + inner(us, vs)*dx
    a2 = lhs(F2)
    L2 = rhs(F2)
    A2 = assemble(a2)
    return A2, L2

def output_data(file_fluid, file_solid, u0, p0, f, disp, force, velocity, t, n):
    print(f"{n}")
    if time_manager.should_output(n):
        print(f"{n}: output...")
        file_fluid.write(u0, t)
        file_fluid.write(p0, t)
        file_fluid.write(f, t)
        file_solid.write(disp, t)
        file_solid.write(force, t)
        file_solid.write(velocity, t)


# Create functions for fluid
u0 =   Function(Vf,   name="velocity")
u0_1 = Function(Vf_1, name="velocity 1st order")
p0 =   Function(Vp,   name="pressure")
f =    Function(Vf_1, name="force")

# Create functions for solid
velocity = Function(Vs, name="velocity")
disp     = Function(Vs, name="displacement")
force    = Function(Vs, name="force")
disp.interpolate(Expression(("x[0]", "x[1]"), degree=2))
ib_interpolation.evaluate_current_points(disp._cpp_object)

# Define fluid solver object
navier_stokes_solver = TaylorHoodSolver(u0, p0, f, dt, nu, stab=stab, alpha=alpha)
bcu, bcp = calculate_fluid_boundary_conditions(navier_stokes_solver.W)

# Define trial and test functions for solid solver
us = TrialFunction(Vs)
vs = TestFunction(Vs)
A2, L2 = calculate_constituitive_model(disp, vs, us)

# Define output path
file_solid_name = unique_filename(os.path.basename(__file__), "note", "/solid.xdmf")
file_fluid_name = unique_filename(os.path.basename(__file__), "note", "/fluid.xdmf")
file_excel_name = unique_filename(os.path.basename(__file__), "note", "/volume.xlsx")
file_log_name   = unique_filename(os.path.basename(__file__), "note", "/info.log")
file_parameters_name   = unique_filename(os.path.basename(__file__), "note", "/parameters.json")
file_solid = create_xdmf_file(solid_mesh.mpi_comm(), file_solid_name)
file_fluid = create_xdmf_file(fluid_mesh.mpi_comm(), file_fluid_name)
logger.add(file_log_name)
logger.info(file_solid_name)
logger.info(file_fluid_name)
logger.info(file_excel_name)
write_paramters(file_parameters_name, beta=1)
print(f"solid_mesh.hmax() {solid_mesh.hmax()}, hmin() {solid_mesh.hmin()}")
print(f"fluid_mesh.hmax() {fluid_mesh.hmax()}, hmin() {fluid_mesh.hmin()}")
print("solid fluid mesh ratio(>2) = ", fluid_mesh.hmin() / solid_mesh.hmax())


t = dt
En = 0.0 # elastic energy
for n in range(1, num_steps+1):
    # step 1. calculate velocity and pressure
    # 计算流体的速度和压力
    navier_stokes_solver.update(u0, p0)
    u1, p1 = navier_stokes_solver.solve(bcu, bcp)
    u0.assign(u1)
    p0.assign(p1)
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
    # output_data(file_fluid, file_solid, u0, p0, f, disp, force, velocity, t, n)
    t = n*dt
    print(t, assemble(inner(u0, u0)*dx))


pe = interpolate(pressure_exact, p0.function_space())
error = calculate_error(pe, p0)

File("pe.pvd") << pe
File("p0.pvd") << p0

logger.info(f"error: {error}")
