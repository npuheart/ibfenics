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
from ibfenics.nssolver import TaylorHoodSolver
from ibfenics.io import unique_filename, create_xdmf_file, write_excel

def advance_disp_bdf2(disp, disp_, velocity, dt):
    Vs = disp.function_space()
    temp_disp = Function(Vs)
    temp_disp.vector()[:] = velocity.vector()[:]*dt*2.0/3.0 + 4.0/3.0*disp.vector()[:] - 1.0/3.0*disp_.vector()[:]
    disp_.vector()[:] = disp.vector()[:]
    disp.vector()[:] = temp_disp.vector()[:]

def advance_disp_be(disp, velocity, dt):
    disp.vector()[:] = velocity.vector()[:]*dt + disp.vector()[:]

# Define time parameters
T = 10.0
dt = 1/8000
num_steps = int(T/dt)

# Define fluid parameters
rho = 1.0
nu = 0.001
n_mesh_fluid = 32

# Define solid parameters
nu_s = 1.0/0.0625 
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

orders       = [order_velocity, order_pressure, order_displacement]
seperations  = [n_mesh_fluid, n_mesh_fluid]
box_points   = [Point(0,0), Point(1, 1)]
circle_outer = Circle(Point(0.5,0.5), 0.25)
circle_inner = Circle(Point(0.5,0.5), 0.25-0.0625)
solid_mesh   = generate_mesh(circle_outer-circle_inner, n_mesh_solid)
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
u0 =   Function(Vf,   name="velocity")
u0_1 = Function(Vf_1, name="velocity 1st order")
p0 =   Function(Vp,   name="pressure")
f =    Function(Vf_1, name="force")

# Create functions for solid
velocity = Function(Vs, name="velocity")
disp     = Function(Vs, name="displacement")
force    = Function(Vs, name="force")

# disp.interpolate(InitialDisplacement())
ib_interpolation.evaluate_current_points(disp._cpp_object)

# Define interpolation object and fluid solver object
navier_stokes_solver = TaylorHoodSolver(u0, p0, f, dt, nu, stab=stab, alpha=alpha)
W = navier_stokes_solver.W

# Define boundary conditions for fluid solver
bcu_1 = DirichletBC(W.sub(0), Constant((0,0)), "near(x[1],1.0)")
bcu_2 = DirichletBC(W.sub(0), Constant((0,0)), "near(x[1],0.0) || near(x[0],0.0) || near(x[0],1.0)")
bcp_1 = DirichletBC(W.sub(1), Constant(0), "near(x[1],0.0) && near(x[0],0.0)", "pointwise")
bcu = [bcu_1, bcu_2]
bcp = [bcp_1]

# Define trial and test functions for solid solver
us = TrialFunction(Vs)
vs = TestFunction(Vs)

# TODO: Define solid constituitive model
# Define variational problem for solid solver
# r = as_vector((-cos(s[0]/R), -sin(s[0]/R)))
# G = mu/omega/(1+s[1])/R*r
# H = - inner(G, V)*dx + inner(U, V)*dx
F = grad(disp)
P = nu_s*(F-inv(F).T)
F2 = inner(P, grad(vs))*dx + inner(us, vs)*dx
a2 = lhs(F2)
L2 = rhs(F2)
A2 = assemble(a2)

# Define output path
if stab:
    note = f"{nu_s}-{dt}-{alpha}"
else :
    note = f"{nu_s}-{dt}"

file_solid_name = unique_filename(os.path.basename(__file__), note, "/solid.xdmf")
file_fluid_name = unique_filename(os.path.basename(__file__), note, "/fluid.xdmf")
file_excel_name = unique_filename(os.path.basename(__file__), note, "/volume.xlsx")
file_log_name   = unique_filename(os.path.basename(__file__), note, "/info.log")
file_solid = create_xdmf_file(solid_mesh.mpi_comm(), file_solid_name)
file_fluid = create_xdmf_file(fluid_mesh.mpi_comm(), file_fluid_name)
logger.add(file_log_name)


# For post-processing
def elastic_energy(disp):
    F = grad(disp)
    C = F*F.T
    Ic = tr(C)
    J = det(F)
    return assemble(0.5*nu_s*(Ic - 2 - 2*ln(J))*dx)


def calculate_volume(X):
    volume_J = assemble(det(grad(X))*dx)
    print("体积：", volume_J)
    return volume_J


t = dt
volume_list = []
En = 0.0 # elastic energy
for n in range(1, num_steps+1):
    En = elastic_energy(disp)+ assemble(0.5*rho*inner(u0, u0)*dx)
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
    if n % 10 == 0:
        file_fluid.write(u0, t)
        file_fluid.write(p0, t)
        file_fluid.write(f, t)
        file_solid.write(disp, t)
        file_solid.write(force, t)
        file_solid.write(velocity, t)
    volume_list.append(calculate_volume(disp))
    t = n*dt
    print(t)

write_excel(volume_list, file_excel_name)

