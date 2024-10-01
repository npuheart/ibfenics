from fenics import *
from mshr import *
from ibfe import Interaction

order_velocity = 2
order_pressure = 1
order_displacement = 1
n_mesh = 60
n_solid_mesh = 60
solid_mesh = generate_mesh(Circle(Point(0.6,0.5), 0.2), n_solid_mesh)
print("Fluid mesh: ", n_mesh, "Solid mesh: ", solid_mesh.num_entities(2), "cells")
orders = [order_velocity, order_pressure, order_displacement]
points = [Point(0,0), Point(1, 1)]
seperations = [n_mesh,n_mesh]
interaction = Interaction(points, seperations, solid_mesh, orders)

ib_interpolation = interaction.ib_interpolation
Vs = interaction.Vs
Vf = interaction.Vf
Vf_1 = interaction.Vf_1
Vp = interaction.Vp


# 测试 Interpolation 和 Distribution
F = interpolate(Expression(("-sin(x[0])*cos(x[1])","cos(x[0])*sin(x[1])"),degree=2), Vs)
f = interpolate(Expression(("-sin(x[0])*cos(x[1])","cos(x[0])*sin(x[1])"),degree=4), Vf_1)

# F1 = interpolate(Expression(("x[0]","x[1]"),degree=2), Vs)
# f1 = interpolate(Expression(("x[0]","x[1]"),degree=4), Vf)

# F2 = interpolate(Expression(("x[0]*x[0]","x[1]*x[1]"),degree=2), Vs)
# f2 = interpolate(Expression(("x[0]*x[0]","x[1]*x[1]"),degree=4), Vf)

# F = F1
# f = f1





EF = Function(Vf_1)
If = Function(Vs)

ib_interpolation.fluid_to_solid(f._cpp_object, If._cpp_object)
ib_interpolation.solid_to_fluid(EF._cpp_object, F._cpp_object)

u = TrialFunction(Vf_1)
v = TestFunction(Vf_1)
a = inner(u,v)*dx
Af = assemble(a)



EF_ = Function(Vf_1)
EF_vector = Function(Vf_1)
weights = []
n = Vf_1.dim()
for i in range(0,n):
    function = Function(Vf_1)
    function.vector()[i] = 1
    weights.append(assemble(function[i%2]*dx))


for i in range(len(weights)):
    # EF_vector.vector()[i] = EF.vector()[i]*weights[i]
    EF_vector.vector()[i] = EF.vector()[i]/n_mesh/n_mesh
# for i in range(len(EF_vector.vector())):
    # EF_vector.vector()[i] = EF.vector()[i]/n_mesh/n_mesh/2/2

solve(Af, EF_.vector(), EF_vector.vector())

u = TrialFunction(Vs)
v = TestFunction(Vs)
a = inner(u,v)*dx
As = assemble(a)

If_ = Function(Vs)
If_vector = Function(Vs)
weights = []
n = Vs.dim()
for i in range(0,n):
    function = Function(Vs)
    function.vector()[i] = 1
    weights.append(assemble(function[i%2]*dx))

for i in range(len(weights)):
    If_vector.vector()[i] = If.vector()[i]*weights[i]

solve(As, If_.vector(), If_vector.vector())

# print(assemble(inner(F, If_)*dx))
print(assemble(inner(EF_, EF_)*dx))
print(assemble(inner(EF, EF)*dx))
print(assemble(inner(EF_, f)*dx))
print(assemble(inner(EF, f)*dx))
print(assemble(inner(EF_, f)*dx) - assemble(inner(F, If_)*dx))
print(assemble(inner(EF, f)*dx)  - assemble(inner(F, If)*dx))

# print(EF(0.5,0.5), EF_(0.5,0.5))
File("EF_.pvd") << EF_
File("EF.pvd") << EF
File("F.pvd") << F
File("f.pvd") << f
File("If_.pvd") << If_
File("If.pvd") << If
print(If(0.5,0.5), If_(0.5,0.5))

import numpy as np
print(np.sqrt(assemble(inner(If_-F, If_-F)*dx)))
print(np.sqrt(assemble(inner(If-F, If-F)*dx)))
# print(np.sqrt(assemble(inner(If-F, If-F)*dx)))

# print(np.sqrt(assemble(div(If)*div(If)*dx)))
# print(np.sqrt(assemble(div(f)*div(f)*dx)))
# print(np.sqrt(assemble(div(If_)*div(If_)*dx)))


# print(np.sqrt(assemble(div(EF)*div(EF)*dx)))
# print(np.sqrt(assemble(div(F)*div(F)*dx)))
# print(np.sqrt(assemble(div(EF_)*div(EF_)*dx)))
# assemble(inner(EF_, u)*dx)

