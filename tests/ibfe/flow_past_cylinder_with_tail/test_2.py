# Copyright (C) 2024 Pengfei Ma
#
# This file is part of ibfenics (https://github.com/npuheart/ibfenics)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
#
# email : ibfenics@pengfeima.cn
#
# brief : flow past cylinder with tail

import os
import numpy as np
from loguru import logger

from mshr import *
from fenics import *

from ibfenics import Interaction, UserIU
from ibfenics.nssolver import TaylorHoodSolver
from ibfenics.io import unique_filename, create_xdmf_file, write_excel


def advance_disp_be(disp, velocity, dt):
    disp.vector()[:] = velocity.vector()[:] * dt + disp.vector()[:]


# Standard units
units = UserIU(g=0.001, cm=0.01)

# Define time parameters
T = 12.0 * units.s
dt = 1 / 20000 * units.s
num_steps = int(T / dt)

# Define fluid parameters
rho = 1.0
nu = 0.001
n_mesh_fluid = 32

# Define solid parameters
n_mesh_solid = 40
beta_s = 1e7
kappa_stab = 1e6
G_s = 1000.0

# Define stablization parameters
alpha = 1.0 * dt
stab = False
delta = 0.1
SAV = 1.0

# Define finite element parameters
order_velocity = 2
order_pressure = 1
order_displacement = 1

# Define spatial discretizations
orders = [order_velocity, order_pressure, order_displacement]
seperations = [n_mesh_fluid * 6, n_mesh_fluid]
box_points = [Point(0, 0), Point(2.2, 0.41)]
from local_mesh import solid_mesh, bdry, domains, dx, marker_circle

interaction = Interaction(box_points, seperations, solid_mesh, orders)

fluid_mesh = interaction.fluid_mesh
ib_mesh = interaction.ib_mesh
ib_interpolation = interaction.ib_interpolation
Vs = interaction.Vs
Vf = interaction.Vf
Vf_1 = interaction.Vf_1
Vp = interaction.Vp

File("b.pvd") << fluid_mesh
print(f"solid_mesh.hmax() {solid_mesh.hmax()}, hmin() {solid_mesh.hmin()}")
print(f"fluid_mesh.hmax() {fluid_mesh.hmax()}, hmin() {fluid_mesh.hmin()}")
print("solid fluid mesh ratio(>2) = ", fluid_mesh.hmin() / solid_mesh.hmax())


# Create functions for fluid
u0 = Function(Vf, name="velocity")
u0_1 = Function(Vf_1, name="velocity 1st order")
p0 = Function(Vp, name="pressure")
f = Function(Vf_1, name="force")

# Create functions for solid
velocity = Function(Vs, name="velocity")
disp = Function(Vs, name="displacement")
force = Function(Vs, name="force")
disp.interpolate(Expression(("x[0]", "x[1]"), degree=2))
ib_interpolation.evaluate_current_points(disp._cpp_object)

# Define interpolation object and fluid solver object
navier_stokes_solver = TaylorHoodSolver(u0, p0, f, dt, nu, stab=stab, alpha=alpha)
W = navier_stokes_solver.W

# Define boundary conditions for fluid solver
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


# Define solid constituitive model
def calculate_constituitive_model(disp, Vs):
    # Define trial and test functions for solid solver
    us = TrialFunction(Vs)
    vs = TestFunction(Vs)
    X0 = SpatialCoordinate(solid_mesh)
    # Define neo-Hookean material
    F = grad(disp)
    I3 = det(F) * det(F)
    P = G_s / 2.0 * (F - inv(F).T) + kappa_stab * ln(I3) * inv(F).T
    F2 = inner(P, grad(vs)) * dx + inner(us, vs) * dx
    # Define stablization term
    F2 += kappa_stab * inner(grad(us), grad(vs)) * dx
    # Define constraints for cylinder
    F2 += beta_s * inner((disp - X0), vs) * dx(subdomain_id=marker_circle)
    a2 = lhs(F2)
    L2 = rhs(F2)
    A2 = assemble(a2)
    return L2, A2


L2, A2 = calculate_constituitive_model(disp, Vs)

# Define output path
note = "note"
file_solid_name = unique_filename(os.path.basename(__file__), note, "/solid.xdmf")
file_fluid_name = unique_filename(os.path.basename(__file__), note, "/fluid.xdmf")
file_excel_name = unique_filename(os.path.basename(__file__), note, "/volume.xlsx")
file_log_name = unique_filename(os.path.basename(__file__), note, "/info.log")
file_solid = create_xdmf_file(solid_mesh.mpi_comm(), file_solid_name)
file_fluid = create_xdmf_file(fluid_mesh.mpi_comm(), file_fluid_name)
logger.add(file_log_name)
print(f"file_solid_name: {file_solid_name}")
print(f"file_fluid_name: {file_fluid_name}")


# For post-processing
def elastic_energy(disp):
    F = grad(disp)
    C = F * F.T
    Ic = tr(C)
    J = det(F)
    return assemble(0.5 * nu_s * (Ic - 2 - 2 * ln(J)) * dx)


def calculate_volume(X):
    volume_J = assemble(det(grad(X)) * dx)
    print("体积：", volume_J)
    return volume_J


t = dt
volume_list = []
En = 0.0  # elastic energy
for n in range(1, num_steps + 1):
    # En = elastic_energy(disp)+ assemble(0.5*rho*inner(u0, u0)*dx)
    # step 1. calculate velocity and pressure
    # 计算流体的速度和压力
    navier_stokes_solver.update(u0, p0)
    u1, p1 = navier_stokes_solver.solve(bcus, bcps)
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
    file_fluid.write(u0, t)
    file_fluid.write(p0, t)
    file_fluid.write(f, t)
    file_solid.write(disp, t)
    file_solid.write(force, t)
    file_solid.write(velocity, t)
    volume_list.append(calculate_volume(disp))
    t = n * dt
    print(t)

write_excel(volume_list, file_excel_name)
