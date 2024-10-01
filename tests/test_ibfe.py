# Copyright (C) 2024 Pengfei Ma
#
# This file is part of ibfenics (https://github.com/npuheart/ibfenics)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
# email : ibfenics@pengfeima.cn

from ibfe import Interaction
from ibfe.nssolver import TaylorHoodSolver
from ibfe.io import unique_filename, create_xdmf_file, write_excel

from dolfin import *
from mshr import *
import os
import numpy as np


# Define time parameters
T = 10.0
dt = 0.001

# Define fluid parameters
rho = 1.0
nu = 0.01
n_mesh_fluid = 32

# Define solid parameters
nu_s = 0.2
n_mesh_solid = 32

# Define stablization parameters
alpha = 1.0*dt
stab = False

# Define finite element parameters
order_velocity = 2 
order_pressure = 1
order_displacement = 1

orders       = [order_velocity, order_pressure, order_displacement]
seperations  = [n_mesh_fluid, n_mesh_fluid]
box_points   = [Point(0,0), Point(1, 1)]
solid_mesh   = generate_mesh(Circle(Point(0.6,0.5), 0.2), n_mesh_solid)
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
disp.interpolate(Expression(("x[0]", "x[1]"), degree=2))

# Create fluid solver
navier_stokes_solver = TaylorHoodSolver(u0, p0, f, dt, nu, stab=stab, alpha=alpha)
W = navier_stokes_solver.W

# Define boundary conditions for fluid solver
noslip = DirichletBC(W.sub(0), (0, 0), "near(x[0],1) || near(x[0],0) || near(x[1],0)")
upflow = DirichletBC(W.sub(0), (1, 0), "near(x[1],1)")
pinpoint = DirichletBC(W.sub(1), 0, "near(x[0],0) && near(x[1],0)", "pointwise")
bcu = [noslip, upflow]
bcp = [pinpoint]

# Define trial and test functions for solid solver
us = TrialFunction(Vs)
vs = TestFunction(Vs)

# Define variational problem for solid solver
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
file_solid = create_xdmf_file(solid_mesh.mpi_comm(), file_solid_name)
file_fluid = create_xdmf_file(fluid_mesh.mpi_comm(), file_fluid_name)


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
num_steps = int(T/dt)
for n in range(1, num_steps+1):
    En = elastic_energy(disp)
    print(En)
    if np.isnan(En):
        print("The simulation is blowing up!")
        break
    print("energy: ", En)
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
    disp.vector()[:] = velocity.vector()[:]*dt + disp.vector()[:]
    ib_interpolation.evaluate_current_points(disp._cpp_object)
    # step 4. calculate body force.
    b2 = assemble(L2)
    solve(A2, force.vector(), b2)
    # step 5. interpolate force from solid to fluid
    ib_interpolation.solid_to_fluid(f._cpp_object, force._cpp_object)
    # step 6. update variables and save to file.
    if n%10 == 0:
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






