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

from ibfenics import Interaction
from ibfenics.nssolver import SAVTaylorHoodSolverBDF2
from ibfenics.io import unique_filename, create_xdmf_file, TimeManager, write_paramters
TaylorHoodSolverBDF2_1 = SAVTaylorHoodSolverBDF2.TaylorHoodSolverBDF2_1
TaylorHoodSolverBDF2_2 = SAVTaylorHoodSolverBDF2.TaylorHoodSolverBDF2_2
modified_energy        = SAVTaylorHoodSolverBDF2.modified_energy
# calculate_SAV_2        = SAVTaylorHoodSolverBDF2.calculate_SAV_2
calculate_SAV          = SAVTaylorHoodSolverBDF2.calculate_SAV

from local_mesh import get_mesh
from ref_coordinates import FiberForce

# Define time parameters
T =  0.00005
dt = 0.00001
num_steps = int(T/dt)
time_manager = TimeManager(T, num_steps, 20)

# Define fluid parameters
nu = 0.01
n_mesh_fluid = 32

# Define solid parameters
n_mesh_solid = 40

# Define stablization parameters
alpha = 1.0*dt
stab  = False
delta = 0.1
SAV   = 1.0

# Define finite element parameters
order_velocity = 2
order_pressure = 1
order_displacement = 1


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

def calculate_fluid_boundary_conditions_sav(W):
    bcp = DirichletBC(W.sub(1), Constant(0), "near(x[1],0.0) && near(x[0],0.0)", "pointwise")
    bcu_3 = DirichletBC(W.sub(0), Constant((0,0)), "on_boundary")
    bcus_2 = [bcu_3]
    bcps_2 = [bcp]
    return bcus_2, bcps_2

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

orders       = [order_velocity, order_pressure, order_displacement]
seperations  = [n_mesh_fluid, n_mesh_fluid]
box_points   = [Point(0,0), Point(1, 1)]
solid_mesh   = get_mesh(n_mesh_solid)
interaction  = Interaction(box_points, seperations, solid_mesh, orders)

fluid_mesh          = interaction.fluid_mesh
ib_mesh             = interaction.ib_mesh
ib_interpolation    = interaction.ib_interpolation
Vs                  = interaction.Vs
Vf                  = interaction.Vf
Vf_1                = interaction.Vf_1
Vp                  = interaction.Vp

print(f"solid_mesh.hmax() {solid_mesh.hmax()}, hmin() {solid_mesh.hmin()}")
print(f"fluid_mesh.hmax() {fluid_mesh.hmax()}, hmin() {fluid_mesh.hmin()}")
print("solid fluid mesh ratio(>2) = ", fluid_mesh.hmin() / solid_mesh.hmax())

# Create functions for fluid
u0_ = Function(Vf, name="velocity_")
p0_ = Function(Vp, name="pressure_")
u0 =   Function(Vf,   name="velocity")
u0_1 = Function(Vf_1, name="velocity 1st order")
p0 =   Function(Vp,   name="pressure")
f =    Function(Vf_1, name="force")

# Create functions for solid
velocity = Function(Vs, name="velocity")
disp     = Function(Vs, name="displacement")
disp_    = Function(Vs, name="displacement_")
force    = Function(Vs, name="force")
disp.interpolate(Expression(("x[0]", "x[1]"), degree=2))
disp_.interpolate(Expression(("x[0]", "x[1]"), degree=2))
ib_interpolation.evaluate_current_points(disp._cpp_object)

# Define interpolation object and fluid solver object
navier_stokes_solver_1 = TaylorHoodSolverBDF2_1(u0, u0, p0, dt, nu)
navier_stokes_solver_2 = TaylorHoodSolverBDF2_2(u0, u0, p0, f, dt, nu)
W = navier_stokes_solver_1.W
bcu_1, bcp_1 = calculate_fluid_boundary_conditions(W)
bcu_2, bcp_2 = calculate_fluid_boundary_conditions_sav(W)

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

t = dt
# volume_list = []
# En = elastic_energy(disp)+ assemble(0.5*rho*inner(u0, u0)*dx)
# qn = np.sqrt(En + delta)
# qnm1 = qn
for n in range(1, num_steps+1):
    # calculate energy and qn
    # En = elastic_energy(disp)+ assemble(0.5*rho*inner(u0, u0)*dx)
    # qnm1 = qn
    # qn = SAV*np.sqrt(En+delta)
    # print(En)
    # if np.isnan(En):
    #     print("The simulation is blowing up!")
    #     break
    # # qn = np.sqrt(modified_energy(u0, En, rho, delta))
    # print("energy: ", En, qn)
    # step 1. calculate velocity and pressure
    # 计算流体的速度和压力，需要计算两个子问题
    navier_stokes_solver_1.update(u0, p0)
    navier_stokes_solver_2.update(u0, p0)
    un   = u0
    unp1 = u0
    unm1 = u0_
    u1, p1 = navier_stokes_solver_1.solve(bcu_1, bcp_1)
    u2, p2 = navier_stokes_solver_2.solve(bcu_2, bcp_2)
    N = FacetNormal(fluid_mesh)
    # SAV = calculate_SAV_2(dt, dt, nu, En, qn, qnm1, unp1, un, unm1, u1, u2, p1, p2, delta, alpha, N, rho)
    # SAV = calculate_SAV(u0, u1, u2, En, dt, nu, delta, qn)
    # print(qn, SAV, n)
    u0_.vector()[:] = u0.vector()[:]
    u0.vector()[:] = u1.vector()[:] + SAV*u2.vector()[:]
    p0.vector()[:] = p1.vector()[:] - SAV*p2.vector()[:]
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
    t = n*dt
    print(t, assemble(inner(u0, u0)*dx))


