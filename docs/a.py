from ibfenics1.cpp import IBMesh3D, IBInterpolation3D
from fenics import *
from mshr import Sphere, generate_mesh





n_mesh = 32
points = [Point(0, 0, 0), Point(1, 1, 1)]
seperations = [n_mesh, n_mesh, n_mesh]
order_velocity = 2
order_displacement = 2



ib_mesh = IBMesh3D(points, seperations, order_velocity)
fluid_mesh = ib_mesh.mesh()
solid_mesh = generate_mesh(Sphere(Point(0.5, 0.5,0.5), 0.2), 30)


Vf = VectorFunctionSpace(fluid_mesh, "P", order_velocity)
Vs = VectorFunctionSpace(solid_mesh, "P", order_displacement)



uf = interpolate(Expression(("x[0]", "x[1]", "x[2]"), degree=2), Vf)
ib_mesh.build_map(uf._cpp_object)

inter = IBInterpolation3D(ib_mesh, solid_mesh)
us = interpolate(Expression(("x[0]", "x[1]", "x[2]"), degree=2), Vs)
inter.evaluate_current_points(us._cpp_object)


# 插值一个点
f = interpolate(Expression(("x[0]*x[0]", "x[1]*x[1]", "x[2]*x[2]"), degree=2), Vf)
a = ib_mesh.evaluate(Point(0.51,0.52,0.53), uf._cpp_object)

# # 插值一个函数
# If = Function(Vs)
# inter.fluid_to_solid(f._cpp_object, If._cpp_object)
# File("us.pvd") << If
# File("uf.pvd") << f


# 分布一个函数


Es = Function(Vf)
us = interpolate(Expression(("x[0]", "x[1]", "x[2]"), degree=2), Vs)

dus = TestFunction(Vs)
Ls = inner(us, dus)*dx
bs = assemble(Ls)
us.vector()[:] = bs[:]

inter.solid_to_fluid(Es._cpp_object, us._cpp_object)
File("uf_.pvd") << Es
File("us_.pvd") << us

