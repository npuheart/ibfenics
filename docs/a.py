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
weights = [1.0 for i in range(0, Vs.dim(), 3)]
inter.evaluate_weights(weights)


Es = Function(Vf)
us = interpolate(Expression(("x[0]", "x[1]", "x[2]"), degree=2), Vs)

dus = TestFunction(Vs)
Ls = inner(us, dus)*dx
bs = assemble(Ls)
us.vector()[:] = bs[:]

inter.solid_to_fluid(Es._cpp_object, us._cpp_object)
File("uf.pvd") << Es
File("us.pvd") << us


# from fenics import *
# from ibfenics1.cpp import IBMesh, IBInterpolation


# class Interaction:
#     def __init__(self, points, seperations, solid_mesh, orders):
#         # 背景网格使用一阶元
#         ib_mesh = IBMesh(points, seperations, 1)
#         fluid_mesh = ib_mesh.mesh()
#         Vf_1 = VectorFunctionSpace(fluid_mesh, "P", 1)
#         # 构造函数空间
#         Vf = VectorFunctionSpace(fluid_mesh, "P", orders[0])
#         Vp = FunctionSpace(fluid_mesh, "P", orders[1])
#         Vs = VectorFunctionSpace(solid_mesh, "P", orders[2])
#         # 初始化插值算子和分布算子
#         uf = interpolate(Expression(("x[0]", "x[1]"), degree=2), Vf_1)
#         us = interpolate(Expression(("x[0]", "x[1]"), degree=2), Vs)
#         ib_mesh.build_map(uf._cpp_object)
#         ib_interpolation = IBInterpolation(ib_mesh, solid_mesh)
#         ib_interpolation.evaluate_current_points(us._cpp_object)
#         n = Vs.dim()
#         weights = []
#         for i in range(0, n, 2):
#             function = Function(Vs)
#             function.vector()[i] = 1
#             weights.append(assemble(function[0] * dx))

#         ib_interpolation.evaluate_weights(weights)
#         # 暴露内部变量
#         self.ib_interpolation = ib_interpolation
#         self.Vf = Vf
#         self.Vf_1 = Vf_1
#         self.Vp = Vp
#         self.Vs = Vs
#         self.ib_mesh = ib_mesh
#         self.fluid_mesh = fluid_mesh
#         self.weights = weights

#         # TODO: Wrappers for Dual Operators
