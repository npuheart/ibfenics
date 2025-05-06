# Copyright (C) 2024 Pengfei Ma
#
# This file is part of ibfenics1 (https://github.com/npuheart/ibfenics1)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
# email : mapengfei@mail.nwpu.edu.cn

from ibfenics1.cpp import IBMesh3D, IBInterpolation3D
from fenics import *

Ne = 32
Nl = 10
order_velocity = 2
order_pressure = 1
order_displacement = 2

ib_mesh = IBMesh3D([Point(0, 0, 0), Point(1, 1, 1)], [Ne, Ne, Ne], order_velocity)
fluid_mesh = ib_mesh.mesh()
solid_mesh = BoxMesh(Point(0.4, 0.4, 0.4), Point(0.6, 0.6, 0.6), Nl, Nl, Nl)

Vf = VectorFunctionSpace(fluid_mesh, "P", order_velocity)
Vs = VectorFunctionSpace(solid_mesh, "P", order_displacement)

# 初始化
uf = interpolate(Expression(("x[0]", "x[1]", "x[2]"), degree=2), Vf)
ib_mesh.build_map(uf._cpp_object)
inter = IBInterpolation3D(ib_mesh, solid_mesh)
us = interpolate(Expression(("x[0]", "x[1]", "x[2]"), degree=2), Vs)
inter.evaluate_current_points(us._cpp_object)


# Interpolate a point (观察到二阶精度)
f = interpolate(Expression(("x[0]*x[0]", "x[1]*x[1]", "x[2]*x[2]"), degree=2), Vf)
a = ib_mesh.evaluate(Point(0.51,0.52,0.53), f._cpp_object)

# Interpolate a function
If = Function(Vs)
inter.fluid_to_solid(f._cpp_object, If._cpp_object)
File("us.pvd") << If
File("uf.pvd") << f


# spread a function
Es = Function(Vf)
us = interpolate(Expression(("x[0]", "x[1]", "x[2]"), degree=2), Vs)

# 1. 先将函数与测试函数做内积
dus = TestFunction(Vs)
Ls = inner(us, dus)*dx
bs = assemble(Ls)
us.vector()[:] = bs[:]

# 2. 然后将内积结果传递给流体
inter.solid_to_fluid(Es._cpp_object, us._cpp_object)
File("uf_.pvd") << Es
File("us_.pvd") << us

