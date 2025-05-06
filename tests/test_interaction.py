# Copyright (C) 2024 Pengfei Ma
#
# This file is part of ibfenics1 (https://github.com/npuheart/ibfenics1)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
# email : ibfenics1@pengfeima.cn


from fenics import *
# from mshr import *
from ibfenics1 import Interaction, Interaction2

order_velocity = 2
order_pressure = 1
order_displacement = 1
n_mesh = 64
n_solid_mesh = 160
solid_mesh = generate_mesh(Circle(Point(0.6, 0.5), 0.2), n_solid_mesh)
print("Fluid mesh: ", n_mesh, "Solid mesh: ", solid_mesh.num_entities(2), "cells")
orders = [order_velocity, order_pressure, order_displacement]
points = [Point(0, 0), Point(1, 1)]
seperations = [n_mesh, n_mesh]
interaction = Interaction2(points, seperations, solid_mesh, orders)

ib_interpolation = interaction.ib_interpolation
Vs = interaction.Vs
Vf = interaction.Vf
Vp = interaction.Vp


# 测试 Interpolation 和 Distribution
F = interpolate(
    Expression(("-sin(x[0])*cos(x[1])", "cos(x[0])*sin(x[1])"), degree=2), Vs
)
f = interpolate(
    Expression(("-sin(x[0])*cos(x[1])", "cos(x[0])*sin(x[1])"), degree=2), Vf
)

EF = Function(Vf)
If = Function(Vs)

ib_interpolation.fluid_to_solid(f._cpp_object, If._cpp_object)
ib_interpolation.solid_to_fluid(EF._cpp_object, F._cpp_object)

File("results/If.pvd") << If
File("results/EF.pvd") << EF
File("results/f.pvd") << f
File("results/F.pvd") << F
