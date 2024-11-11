from fenics import *
from ibfenics1.cpp import IBMesh, IBInterpolation


class Interaction:
    def __init__(self, points, seperations, solid_mesh, orders):
        # 背景网格使用一阶元
        ib_mesh = IBMesh(points, seperations, 1)
        fluid_mesh = ib_mesh.mesh()
        Vf_1 = VectorFunctionSpace(fluid_mesh, "P", 1)
        # 构造函数空间
        Vf = VectorFunctionSpace(fluid_mesh, "P", orders[0])
        Vp = FunctionSpace(fluid_mesh, "P", orders[1])
        Vs = VectorFunctionSpace(solid_mesh, "P", orders[2])
        # 初始化插值算子和分布算子
        uf = interpolate(Expression(("x[0]", "x[1]"), degree=2), Vf_1)
        us = interpolate(Expression(("x[0]", "x[1]"), degree=2), Vs)
        ib_mesh.build_map(uf._cpp_object)
        ib_interpolation = IBInterpolation(ib_mesh, solid_mesh)
        ib_interpolation.evaluate_current_points(us._cpp_object)
        n = Vs.dim()
        weights = []
        for i in range(0, n, 2):
            function = Function(Vs)
            function.vector()[i] = 1
            weights.append(assemble(function[0] * dx))

        ib_interpolation.evaluate_weights(weights)
        # 暴露内部变量
        self.ib_interpolation = ib_interpolation
        self.Vf = Vf
        self.Vf_1 = Vf_1
        self.Vp = Vp
        self.Vs = Vs
        self.ib_mesh = ib_mesh
        self.fluid_mesh = fluid_mesh
        self.weights = weights

        # TODO: Wrappers for Dual Operators


class Interaction2:
    def __init__(self, points, seperations, solid_mesh, orders):
        # 背景网格使用一阶元
        order_velocity, order_pressure, order_displacement = orders
        ib_mesh = IBMesh(points, seperations, order_velocity)
        fluid_mesh = ib_mesh.mesh()
        # 构造函数空间
        Vf = VectorFunctionSpace(fluid_mesh, "P", order_velocity)
        Vp = FunctionSpace(fluid_mesh, "P", order_pressure)
        Vs = VectorFunctionSpace(solid_mesh, "P", order_displacement)
        # 初始化插值算子和分布算子
        uf = interpolate(Expression(("x[0]", "x[1]"), degree=2), Vf)
        us = interpolate(Expression(("x[0]", "x[1]"), degree=2), Vs)
        ib_mesh.build_map(uf._cpp_object)
        ib_interpolation = IBInterpolation(ib_mesh, solid_mesh)
        ib_interpolation.evaluate_current_points(us._cpp_object)
        n = Vs.dim()
        weights = []
        for i in range(0, n, 2):
            function = Function(Vs)
            function.vector()[i] = 1
            weights.append(assemble(function[0] * dx))

        ib_interpolation.evaluate_weights(weights)
        # 暴露内部变量
        self.ib_interpolation = ib_interpolation
        self.Vf = Vf
        self.Vp = Vp
        self.Vs = Vs
        self.ib_mesh = ib_mesh
        self.fluid_mesh = fluid_mesh
        self.weights = weights

        # TODO: Wrappers for Dual Operators
