# Copyright (C) 2024 Pengfei Ma
#
# This file is part of ibfenics (https://github.com/ibfenics)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

from ibfe import Interaction
from ibfe.nssolver import TaylorHoodSolver

from dolfin import *
from mshr import *
from datetime import datetime
import os
import numpy as np


# Define fluid solver
rho = 1.0
T = 10.0
n_mesh = 80
dt = 0.001
nu = 0.01
nu_s = 0.2

# stabliization parameter
alpha = 1.0*dt
stab = False

current_file_name = os.path.basename(__file__)
note =  os.path.splitext(current_file_name)[0]
file_id = "data/" + note + "-" + datetime.now().strftime('%Y%m%d-%H%M%S')

order_velocity = 2 
order_pressure = 1
order_displacement = 1
solid_mesh = generate_mesh(Circle(Point(0.6,0.5), 0.2), n_mesh)
orders = [order_velocity, order_pressure, order_displacement]
points = [Point(0,0), Point(1, 1)]
seperations = [64, 64]
interaction = Interaction(points, seperations, solid_mesh, orders)


fluid_mesh = interaction.fluid_mesh
ib_mesh = interaction.ib_mesh
ib_interpolation = interaction.ib_interpolation
Vs = interaction.Vs
Vf = interaction.Vf
Vf_1 = interaction.Vf_1
Vp = interaction.Vp


print("solid_mesh.hmax() ", solid_mesh.hmax())
print("solid_mesh.hmin() ", solid_mesh.hmin())
print("fluid_mesh.hmax() ", fluid_mesh.hmax())
print("fluid_mesh.hmin() ", fluid_mesh.hmin())

print("ratio(>2) = ", fluid_mesh.hmin() / solid_mesh.hmax())

# Define trial and test functions for solid
us = TrialFunction(Vs)
vs = TestFunction(Vs)

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

# Define interpolation object and fluid solver object
navier_stokes_solver = TaylorHoodSolver(u0, p0, f, dt, nu, stab=stab, alpha=alpha)
W = navier_stokes_solver.W

# Define boundary conditions
noslip = DirichletBC(W.sub(0), (0, 0), "near(x[0],1) || near(x[0],0) || near(x[1],0)")
upflow = DirichletBC(W.sub(0), (1, 0), "near(x[1],1)")
pinpoint = DirichletBC(W.sub(1), 0, "near(x[0],0) && near(x[1],0)", "pointwise")
bcu = [noslip, upflow]
bcp = [pinpoint]

F = grad(disp)
P = nu_s*(F-inv(F).T)
# Define variational problem for solid
F2 = inner(P, grad(vs))*dx + inner(us, vs)*dx
a2 = lhs(F2)
L2 = rhs(F2)
A2 = assemble(a2)

# Output Directory name
directory = file_id + "-" + str(nu_s) + "-" + str(dt)
if stab:
    directory = directory + "-" + str(alpha)

print("output directory : ", directory)

# Create files for storing solution
file_solid = XDMFFile(directory+"/solid.xdmf")
file_solid.parameters['rewrite_function_mesh'] = False
file_solid.parameters["functions_share_mesh"] = True
file_solid.parameters["flush_output"] = True

file_fluid = XDMFFile(directory+"/fluid.xdmf")
file_fluid.parameters['rewrite_function_mesh'] = False
file_fluid.parameters["functions_share_mesh"] = True
file_fluid.parameters["flush_output"] = True

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

import pandas as pd

print(volume_list)
df = pd.DataFrame(volume_list)

# 指定要保存的文件名和表单名称
excel_file1 = 'volume_64-fe.xlsx'
sheet_name = 'v'

# 将数据写入 Excel 文件
df.to_excel(excel_file1, sheet_name=sheet_name, index=False)

